#-----------------------------------------------------------------------
# CMAES_Model
#-----------------------------------------------------------------------
#   - center::Vector      Model center
#   - C::Matrix           Covariance matrix
#   - s::Vector           "center" path (computed using shaped noise only)
#   - s_σ::Real           stepsize "path" (single value state about stepsize updated and carried over from one generation to the next)
#   - σ::Real             stepsize
#   - parms::Model_Parms  constant parameters used by the model for updating or setup
#-----------------------------------------------------------------------
# uses  Model_Parms
# N::Integer        chromosome length
# μ::Integer        population size after truncation selection
# λ::Integer        original "sampling" populaiton size
# μ_eff::Real
# c_c::Real
# c_σ::Real
# c_1::Real         one-update proportion
# c_μ::Real         μ-update proportion
# d_σ::Real
# chi_mean::Real    expected value of the chi distribution (expected value of the L2 norm applied to a normal distribution)
# w::Weights        weights used for weighted averages across the population based on the members fitness
#-----------------------------------------------------------------------

function Model_Parms(n::Integer, f::RealFitness;
                     λ = 4 + floor(Int, 3 * log(n)),
                     μ = floor(Int, λ/2),
                     μ_ = μ + 0.5,            # originally (λ + 1)/2 -> changed because we might not want μ = λ/2
                     w = normalize(map((i)->(log(μ_) - log(i)), 1:μ), 1),
                     #w = fill(1/μ,μ),
                     #μ_eff = 1 / sum(map(i->(abs(w[i])^(2-(log(μ)/n))),1:μ)),
                     μ_eff = 1 / sum(map((i)->(w[i])^2,1:μ)),
                     c_c = (4 + μ_eff/n) / (n + 4 + 2μ_eff/n),
                     c_σ = (μ_eff + 2) / (n + μ_eff + 5),
                     c_1 = 2 / ((n + 1.3)^2 + μ_eff),
                     α_μ = 2,
                     c_μ = min(1 - c_1, α_μ * (μ_eff - 2 + 1/μ_eff) / ((n + 2)^2 + α_μ * μ_eff / 2)),
                     d_σ = 1 + 2 * max(0, sqrt((μ_eff - 1) / (n + 1))) + c_σ,
                     #d_σ = (0.7 * (1/c_σ)) + 1,
                     chi_mean = sqrt(n) * (1 - 1/4n + 1/21n^2),
                     orig_scale = 1.0,
                     best_scale = 1.0)
  Model_Parms(n, λ, μ, μ_eff, c_c, c_σ, c_1, c_μ, d_σ, chi_mean, Weights(w), direction(f), orig_scale, best_scale)
end

#-------------------------
#  Model internal update functions
#-------------------------

#---------------
# Internal Update helper functions
function lambda_scale_parms(m::CMAES_Model)
    (m.parms.orig_scale, m.parms.best_scale)
end

function covarmatrix_updateparms(m::CMAES_Model)
  (m.C, m.p_c, m.h_σ, m.parms.c_c, m.parms.c_1, m.parms.c_μ, m.parms.μ_eff)
end

function stepsize_updateparms(m::CMAES_Model)
  (m.σ, m.p_σ, m.parms.c_σ, m.parms.d_σ, m.parms.μ_eff, m.parms.chi_mean)
end

#=function stepsize_updateparms_(m::CMAES_Model)
  (m.σ_, m.p_σ_, m.parms.c_σ, m.parms.d_σ, m.parms.μ_eff, m.parms.chi_mean)
end=#

function setcovarmatrix!(model::CMAES_Model, p_c, C)
  model.p_c = p_c
  model.C = C
end

function setstepsize!(model::CMAES_Model, p_σ, σ)
  model.p_σ = p_σ
  model.σ = σ
end

#=function setstepsize!_(model::CMAES_Model, p_σ_, σ_)
  model.p_σ_ = p_σ_
  model.σ_ = σ_
end=#
#---------------
# Internal Update functions

#  forms recombinant single intermediate parent for next generation
function center!(m::CMAES_Model, w)
    if any(i -> isnan(i), w)
        m.center = fill(NaN,size(m.center,1))
    else
        m.center += m.σ * w
    end
end

#update center_ with best solution from sorted shadow population
center!_(state::CMAES_State,v::Vector{Float64}) = state.nModel_shadow.center_ = v


function h_σ(p_σ, gen::Integer, parms::Model_Parms)
  (c_σ, chi_mean, n) = (parms.c_σ, parms.chi_mean, parms.N)
  if any(i -> isnan(i), p_σ)
      NaN
  else
      (norm(p_σ) / sqrt(1 - (1 - c_σ)^(2(gen + 1))) < (1.4 + 2/(n + 1)) * chi_mean) ? 1.0 : 0.0
  end
end

h_σ(model::CMAES_Model) = model.h_σ

function h_σ!(m::CMAES_Model, gen::Integer)
  m.h_σ = h_σ(m.p_σ, gen, m.parms)
end

function covarmatrix!(model::CMAES_Model, W::Noise)
  (C, p_c, h_σ, c_c, c_1, c_μ, μ_eff) = covarmatrix_updateparms(model)
  if any(i -> isnan(i), members(W))
      p_c = fill(NaN, size(p_c,1))
      C = fill(NaN, size(C))
  else
      p_c = (1 - c_c) * p_c + h_σ * sqrt(c_c * (2 - c_c) * μ_eff) * weightedavg(W)            # update path
      C = (1 - c_1 - c_μ) * C + c_1 * p_c * p_c' + c_μ * covariance(W, weights(model))        # update covariance matrix
  end
  setcovarmatrix!(model, p_c, C)
end

function stepsize!(model::CMAES_Model, w)
  (σ, p_σ, c_σ, d_σ, μ_eff, chi_mean) = stepsize_updateparms(model)
  if any(i -> isnan(i),w)
      p_σ = fill(NaN, size(p_σ,1))
      σ = NaN
  else
      p_σ = (1 - c_σ) * p_σ + sqrt(c_σ * (2 - c_σ) * μ_eff) * invsqrtC(model) * w  #stepsize path update
      σ = σ * exp(c_σ * (norm(p_σ)/chi_mean - 1) / d_σ) # update stepsize
  end
  setstepsize!(model, p_σ, σ)
end

#=function updateOtherStepSize!(model::CMAES_Model, w)
    (σ_, p_σ_, c_σ, d_σ, μ_eff, chi_mean) = stepsize_updateparms_(model)
    p_σ_ = (1 - c_σ) * p_σ_ + sqrt(c_σ * (2 - c_σ) * μ_eff) * invsqrtC(model) * w
    σ_ = σ_ * exp(c_σ * (norm(p_σ_)/chi_mean - 1)/ d_σ)
    setstepsize!_(model, p_σ_, σ_)
end=#


# eig(m.C) -> (v, B) where v are the eigenvalues as a vector and B are the eigenvectors
# if we set V = Diag(v) then BVB' = C
# if we set D = sqrt(V) we get BDDB' = C
# if (sum(D_sq .< 0) > 0) println("negative eigenvalue: eigenvalues = $D_sq") end

eigendecomp(m::CMAES_Model) = (m.γ, m.B)
sqrtC(m::CMAES_Model) = m.B * m.D * m.B'
invsqrtC(m::CMAES_Model) = m.B * m.I * m.B'

zeroeigenval(m::CMAES_Model) = any(m.γ .== 0.0)

function invsqrtC!(m::CMAES_Model)
  if zeroeigenval(m)
    throw(ZeroEigError())
  else
    m.I = inv(m.D)
  end
end

function eigendecomp!(m::CMAES_Model)
  try
    (m.γ, m.B) = eigen(m.C)
  catch
    throw(ComplexEigError())
  end
end

negeigenval(m::CMAES_Model) = any(m.γ .< 0.0)

function sqrtC!(m::CMAES_Model)
  if negeigenval(m)
    throw(NegEigError())
  else
    m.D = Diagonal(sqrt.(m.γ))
  end
end

#take best solution as the center, w[1] = 1.0, w[2:end] = 0.0
function w_1!(model::CMAES_Model)
    model.parms.w = Weights(fill(0.0, mu(model)))
    model.parms.w[1] = 1.0
end

c_c!(model::CMAES_Model) = model.parms.c_c = (4 + model.parms.μ_eff/N(model)) / (N(model) + 8 + 2 * model.parms.μ_eff/N(model))
#c_σ!(model::CMAES_Model) = model.parms.c_σ = (model.parms.μ_eff + 2) / (N(model) + model.parms.μ_eff + 5)
function getindex(m::CMAES_Model, choice::Symbol, i)
  if choice == :eigvec
    m.B[:, i]
  elseif choice == :eigval
    m.γ[i]
  else
    error("choice $choice not valid - only :eigvec or :eigval used")
  end
end

#-------------------------
#  Model public functions
#-------------------------

N(model::CMAES_Model) = model.parms.N
center(model::CMAES_Model) = model.center
center_(model::CMAES_Model) = model.center_
center(model::CMAES_Model, popnSize) =
          RegularPopulation(center(model), popnSize; direction = model.parms.direction)
centermember(model::CMAES_Model, objfun::Function) =
          Member(center(model), objfun)
centermember_(model::CMAES_Model, objfun::Function) =
          Member(center_(model), objfun)
σ_estimate(model::CMAES_Model) = model.σ
sigma(model::CMAES_Model) = model.σ
#sigma_(model::CMAES_Model) = model.σ_
#sigma!_(model::CMAES_Model, σ_::Float64) = model.σ_ = σ_
covar(model::CMAES_Model) = model.C
+(m::CMAES_Model, n::Noise) = n + m
weights(model::CMAES_Model) = model.parms.w
direction(model::CMAES_Model) = model.parms.direction
mu(model::CMAES_Model) = model.parms.μ
lambda(model::CMAES_Model) = model.parms.λ
chrlength(model::CMAES_Model) = model.parms.N
μ_eff(model::CMAES_Model) = model.parms.μ_eff
c_σ(model::CMAES_Model) = model.parms.c_σ
c_c(model::CMAES_Model) = model.parms.c_c
c_μ(model::CMAES_Model) = model.parms.c_μ
c_1(model::CMAES_Model) = model.parms.c_1
d_σ(model::CMAES_Model) = model.parms.d_σ
p_σ(model::CMAES_Model) = model.p_σ
p_c(model::CMAES_Model) = model.p_c
α_μ(model::CMAES_Model) = model.parms.α_μ
chi_mean(model::CMAES_Model) = model.parms.chi_mean
orig_scale!(model::CMAES_Model, f::Float64) = model.parms.orig_scale = f
best_scale!(model::CMAES_Model, f::Float64) = model.parms.best_scale = f
orig_scale(model::CMAES_Model) = model.parms.orig_scale
best_scale(model::CMAES_Model) = model.parms.best_scale

#breaking up results of equations
#-------------SIGMA STEP SIZE PATH UPDATE--------------#
#(1-c_σ)p_σ
p_σ_part1(model::CMAES_Model) = (1 - c_σ(model)) * (p_σ(model))
#sqrt(c_σ(2 - c_σ)μ_eff)
p_σ_part2(model::CMAES_Model) = sqrt(c_σ(model)*(2-c_σ(model))*μ_eff(model))
#sqrt(c_σ(2 - c_σ)μ_eff)BD^-1B^T⟨y⟩
p_σ_part3(model::CMAES_Model,n::ShapedNoise) = p_σ_part2(model) * (invsqrtC(model) * weightedavg(n))
#------SIGMA UPDATE-----------#
#c_σ/d_σ
σ_part1(model::CMAES_Model) = (c_σ(model) / d_σ(model))
#||p_σ|| a.k.a norm(p_σ)
σ_part2(model::CMAES_Model) = norm(p_σ(model))
#E||N(0,I)|| a.k.a chi_mean approximation
σ_part3(model::CMAES_Model) = chi_mean(model)
#(norm(p_σ)/norm(chi_mean)) - 1
σ_part4(model::CMAES_Model) = (σ_part2(model) / σ_part3(model)) - 1
#(c_σ/d_σ((norm(p_σ)/norm(chi_mean)) - 1))
σ_part5(model::CMAES_Model) = (σ_part1(model) * σ_part4(model))
#exp((c_σ/d_σ((norm(p_σ)/norm(chi_mean)) - 1)))
σ_part6(model::CMAES_Model) = exp(σ_part5(model))
#-------COVARIANCE PATH UPDATE----------#
#(1-c_c)p_c
p_c_part1(model::CMAES_Model) = (1 - c_c(model)) * p_c(model)
#h_σ*sqrt(c_c(2-c_c)μ_eff)
p_c_part2(model::CMAES_Model) = h_σ(model) * sqrt(c_c(model)*(2-c_c(model))*μ_eff(model))
#h_σ*sqrt(c_c(2-c_c)μ_eff)⟨y⟩
p_c_part3(model::CMAES_Model,n::ShapedNoise) = p_c_part2(model) * weightedavg(n)
#-------COVARIANCE MATRIX UPDATE--------#
#(1 - c_1 - c_μ)C
C_part1(model::CMAES_Model) = (1 - c_1(model) - c_μ(model)) * covar(model)
#c_1(p_c * p_c')
C_part2(model::CMAES_Model) = c_1(model) * (p_c(model) * p_c(model)')
# mu-update - covariance(W, weights(model))
C_part3(model::CMAES_Model, n::Noise) = c_μ(model) * covariance(n, weights(model))

# note: order of the updates is important. Each update is dependant on the previous updates having been run.
function update!(model::CMAES_Model, W::ShapedNoise, gen::Integer)
    center!(model, weightedavg(W))
    stepsize!(model, weightedavg(W))
    h_σ!(model, gen)
    #c_c!(model)
    covarmatrix!(model, W)
end

function update(model::CMAES_Model, W::ShapedNoise, gen::Integer)
    model = deepcopy(model)
    update!(model, W, gen)
    model
end

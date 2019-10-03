#-----------------------------------------------------------------------
#  abstract Noise
#-----------------------------------------------------------------------

# SphericalNoise and ShapedNoise Inner Constructor helper functions

function Noise(n::Noise, members)
  n.values = RegularPopulation(members)
  n.flattened = false
  n
end

function Noise(n::Noise, popn::Population)
  n.values = popn
  n.flattened = false
  n
end

function Noise(n::Noise, members, structure::SortStructure)
  n.values = SortedPopulation(members, structure)
  n.flattened = false
  n
end

function Noise(n::Noise, popn::Population, structure::SortStructure)
  n.values = SortedPopulation(popn, structure)
  n.flattened = false
  n
end

function Noise(n::Noise, members, structure::SortStructure, model::CMAES_Model)
  n.values = SortedPopulation(members, structure)
  flatten!(n, weights(model))
  n.flattened = true
  n
end

function Noise(n::Noise, popn::Population, structure::SortStructure, model::CMAES_Model)
  n.values = SortedPopulation(popn, structure)
  flatten!(n, weights(model))
  n.flattened = true
  n
end

RegularPopulation(n::Noise, members) = RegularPopulation(n.values, members)
RegularPopulation(n::Noise, members, fitness) = RegularPopulation(n.values, members, fitness)

# -------------------------------

# CenterNoise and ZeroNoise acts like constructors
#    - creates noise i.e. error of the center value, i.e. 0
#    - CenterNoise produces a population with a single member with 0's for all loci
#    - ZeroNoise produces a population with λ members, each with 0's for all loci
#    - converts that to a Noise object of the same type that was passed in.
zeronoise(N, λ) = zeros(N, λ)
centernoise(N) = zeros(N, 1)
ZeroNoise(noiseType::Type, N, λ) = noiseType(zeronoise(N, λ))
CenterNoise(noiseType::Type, N) = ZeroNoise(noiseType, N, 1)
ZeroNoise(noise::Noise, λ) = typeof(noise)(zeronoise(chrlength(noise), λ))
CenterNoise(noise::Noise) = ZeroNoise(noise, 1)


# -------------------------------
# Inheritable Functions (can be overridden)

noisevalues(n::Noise) = n.values
members(n::Noise) = members(noisevalues(n::Noise))
popnsize(n::Noise) = popnsize(noisevalues(n))
sorted(n::Noise) = n.sortedValues
issorted(n::Noise) = (typeof(n.sortedValues) == SortedPopulation)
flattened(n::Noise) = n.weightedAvg
isflattened(n::Noise) = n.flattened
weightedavg(n::Noise) = n.weightedAvg
popnsize(n::Noise) = popnsize(n.values)
chrlength(n::Noise) = chrlength(n.values)
getindex(n::Noise, index) = n.values[index]
getindex(n::Noise, choice::Symbol, index) = n.values[choice, index]
setindex!(n::Noise, value, index) = (n.values[index] = value)
setindex!(n::Noise, value, choice::Symbol, index) = (n.values[choice, index] = value)
pcat(n1::Noise, n2::Noise) = typeof(n1)(mcat(noisevalues(n1), noisevalues(n2)))
pcat(n1::Noise, n2::Noise, n3::Noise) = typeof(n1)(mcat(noisevalues(n1), noisevalues(n2), noisevalues(n3)))
+(n::Noise, popn::Population) = noisevalues(n) + popn
+(popn::Population, n::Noise) = popn + noisevalues(n)
+(n1::Noise, n2::Noise) = typeof(n1)(mcat(noisevalues(n1) , noisevalues(n2)))

# important - adds noise to a model to create a regular population
+(n::Noise, m::CMAES_Model) = RegularPopulation(n, σ_estimate(m) * members(n) .+ center(m))

function sort!(n::Noise, structure::SortStructure)
  n.values = SortedPopulation(n.values, structure)
end

function sort(n::Noise, structure::SortStructure)
  sn = deepcopy(n)
  sort!(sn, structure)
  sn
end

function matchedinterleave(n1::Noise, n2::Noise, slotSize, pOrder, cOrder)
  values = matchedinterleave(n1.values, n2.values, slotSize, pOrder, cOrder)
  typeof(n1)(values)
end

function matchedtruncation(n::Noise, indexOrder::Array)
  values = matchedtruncation(n.values, indexOrder)
  typeof(n)(values)
end

function flatten!(n::Noise, w::Weights)
  n.flattened = true
  #if any(i -> isnan(i),n[:chr,:])
    #println("found -Infs in noise")
    #n.weightedAvg = fill(NaN,size(n[:chr,:],1))
   n.weightedAvg = vec(mean(n[:chr, :], w, dims=2))
  #end
end

#=function flatten_with_σ!(n::Noise, w::Weights, model::CMAES_Model)
    n.flattened = true
    size_ = lambda(model) # λ
    size_1 = floor(Int64, size_/2)
    size_2 = size_ - size_1
    #println("μ = $(mu(model)), λ = $(lambda(model))")
    #println("size = $(size_), size_1 = $(size_1), size_2 = $(size_2)")
    orig_λ = floor(Int64, orig_scale(model) * size_1)
    best_λ = size_ - orig_λ  #leftover lambda
    noise = zeros(N(model), mu(model)) # N x λ
    #println("orig_λ = $(orig_λ), best_λ = $(best_λ)")
    #println("length of index = $(length(noisevalues(n).index))")
    for i = 1:mu(model)   # 1:μ
        if  noisevalues(n).index[i] > orig_λ
            noise[:,i] = model.σ_ * n[:chr, i]
        else
            noise[:,i] = model.σ * n[:chr, i]
        end
    end
    #weightedAvg already multiplied by the corresponding sigmas
    #so its all ready to added to the center for center! update without multiplying it by sigma
    #println("dims of weightedAvg = $(size(noise,1)) x $(size(noise,2))")
    n.weightedAvg = vec(mean(noise[:,:], w, dims=2))
    #println("σ_ = $(sigma_(model)), weightedAvg after sigmas = $(n.weightedAvg)")
end=#

function flatten(noise::Noise, w::Weights)
  n = deepcopy(noise)
  flatten!(n, w)
  n
end

function update!(n::Noise, structure::SortStructure, fitnessWeights::Weights)
  sort!(n, structure)
  flatten!(n, fitnessWeights)
end

function covariance(noise::Noise, w::Weights)
  n = chrlength(noise.values)
  covar = zeros(n, n)
  for i = 1:popnsize(noise)
    covar = covar + w[i] * noise[:chr, i] * noise[:chr, i]'
  end
  covar
end

# -------------------------------
#  Conversion functions between two main types of Gaussian noise - spherical and shaped

shapedtosphere(shaped, m::CMAES_Model) = (inv(m.D) * inv(m.B) * shaped)
spheretoshaped(sphere, m::CMAES_Model) = (m.B * m.D * sphere)

#-----------------------------------------------------------------------
#  SphericalNoise <: Noise   # aka E
#-----------------------------------------------------------------------
#    values::Population
#    weightedAvg::Vector
#    flattened::Bool
#-----------------------------------------------------------------------
# Internal constructors
#   - SphericalNoise(members)
#   - SphericalNoise(members, structure::SortStructure)
#   - SphericalNoise(members, structure::SortStructure, model::CMAES_Model)
#-----------------------------------------------------------------------

# -------------------------------
# constructors (and converter from ShapedNoise)

function SphericalNoise(N::Integer, λ::Integer)
  randomMembers = randn(N, λ)
  SphericalNoise(randomMembers)
end

function SphericalNoise(shaped::ShapedNoise)
  ShapedNoise(members(shaped))
end

function SphericalNoise(shaped::ShapedNoise, model::CMAES_Model)
  # needs to be well tested
  sphere = SphericalNoise(shaped)
  for i = 1:popnsize(sphere)
    sphere[:chr, i] = shapedtosphere(shaped[:chr, i], model)
  end
  sphere
end

# -------------------------------
# Public functions

# same as Noise - no overides


#-----------------------------------------------------------------------
#  ShapedNoise <: Noise   # aka W
#-----------------------------------------------------------------------
#    *** same instance variables as SphericalNoise ***
#    values::Population
#    weightedAvg::Vector
#    flattened::Bool
#-----------------------------------------------------------------------
# Internal constructors
#   - ShapedNoise(members)
#   - ShapedNoise(members, structure::SortStructure)
#   - ShapedNoise(members, structure::SortStructure, model::CMAES_Model)
#-----------------------------------------------------------------------

# -------------------------------
# constructors (and converter from SphericalNoise)

function ShapedNoise(sphere::SphericalNoise)
  ShapedNoise(members(sphere))
end

function ShapedNoise(sphere::SphericalNoise, model::CMAES_Model; dualcenter = false)
    if !dualcenter
        shaped = ShapedNoise(sphere)
        for i=1:popnsize(shaped)
            shaped[i] = spheretoshaped(sphere[:chr,i], model)
        end
        return shaped
    else

        shaped = ShapedNoise(sphere) # λ
        #println("μ = $(mu(model)), λ = $(lambda(model))")
        size = popnsize(shaped) # λ
        size_1 = floor(Int64, size/2)
        size_2 = size - size_1
        #println("size = $(size)")
        #println("size_2 = $(size_2)")
        if orig_scale(model) == 2.0
            orig_λ = size
            best_λ = 0
        else
            orig_λ = floor(Int64, orig_scale(model) * size_1)
            best_λ = size - orig_λ  #leftover lambda
        end
        shaped_orig_container = deepcopy(shaped)
        shaped_best_container = deepcopy(shaped)

        #println("orig_λ = $(orig_λ), best_λ = $(best_λ)")
        if orig_λ > 0
            shaped_orig = ShapedNoise(members(shaped[:chr,1:orig_λ]))
            for i = 1:orig_λ
                shaped_orig[i] = spheretoshaped(sphere[:chr, i], model)
            end
            shaped_orig_container = shaped_orig
            if orig_λ == size  shaped_best_container = NaN end
        end

        if best_λ > 0
            shaped_best = ShapedNoise(members(shaped[:chr,orig_λ+1:size]))
            a = 0
            for i = (orig_λ+1):size
                a += 1
                shaped_best[a] = spheretoshaped(sphere[:chr, i], model)
            end
            shaped_best_container = shaped_best
            if best_λ == size  shaped_orig_container = NaN end
        end
        if orig_λ == size  shaped_best_container = NaN end
        if best_λ == size  shaped_orig_container = NaN end
        #println("size of orig population = $(popnsize(shaped_orig_container)) & best population = $(popnsize(shaped_best_container))")
        (shaped_orig_container, shaped_best_container)
    end
end
# is σ_estimate accurate enough because it may have hardly changed from the previous generation?
ShapedNoise(popn::Population, model::CMAES_Model; shadow = false)  = ((shadow) ? ShapedNoise((members(popn) .- center(model)) / σ_estimate(model))
                                                                            : ShapedNoise((members(popn) .- center(model)) / σ_estimate(model)))
#ShapedNoise(popn::Population, model::CMAES_Model) = ShapedNoise(members(popn) .- center(model) / σ_estimate(model)) <-shadow


# -------------------------------
# Public functions

# same as Noise - no overides

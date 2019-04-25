#-------------------------
#  CMAES Types
#-------------------------
#  - CMAES_Model
#  - Reproduction_Parms aka Model_Parms
#  - Noise
#  - SphericalNoise, ShapedNoise <: Noise
#  - Selection_Parms
#  - CMA_ES_System, IPop_System         <: CMAES_System     <: System
#  - CMA_ES_State, IPop_State           <: CMAES_State      <: State
#  - CMA_ES_ReturnInfo, iPop_ReturnInfo <: CMAES_ReturnInfo <: EA_ReturnInfo
#  - RealFitness = UnBoundedFitness{Float64}
#-------------------------

const RealFitness = UnBoundedFitness{Float64}


#-------------------------
#  CMAES_Model, including Reproduction_Parms aka Model_Parms

mutable struct Model_Parms
  N::Integer
  λ::Integer
  μ::Integer
  μ_eff::Real
  c_c::Real
  c_σ::Real
  c_1::Real
  c_μ::Real
  d_σ::Real
  chi_mean::Real
  w::Weights
  direction::Symbol
  orig_scale::Float64 # used to scale how much solutions are generated off of center and center_
  best_scale::Float64 # used to scale how much solutions are generated off of center and center_
end

const Reproduction_Parms = Model_Parms

struct ZeroEigError <: Exception end
struct NegEigError <: Exception end
struct ComplexEigError <: Exception end

mutable struct CMAES_Model
  center::Vector{Float64}     # center of the model / population / multi-variate normal distribution
  center_::Vector{Float64}    # center 2 of the model which always gets replaced with the best solution of the previous population
  C::Array{Float64, 2}        # covariance matrix
  B::Array{Float64, 2}        # normalized eigenvectors of C (a matrix)
  D::Diagonal{Float64}        # sq roots of the eigenvalues diagonalized (a matrix)
  I::Diagonal{Float64}        # inv sq roots of the eigenvalues diagonalized (a matrix)
  γ::Vector{Float64}          # eigenvalues of C (a vector)
  p_c::Vector{Float64}        # covariance path
  p_σ::Vector{Float64}        # step-size path
  h_σ::Float64                # used as a switch (makes sure p_c not too big) in the update equations (1.0 if true 0.0 if false)
  σ::Float64                  # estimated standard deviation - aka step-size
  parms::Model_Parms

  function CMAES_Model(parms::Model_Parms, center_init, σ_init; gen = 0)
    N = parms.N
    m = new()
    m.parms = parms
    m.center = center_init
    m.center_ = deepcopy(center_init)
    m.C = Matrix(1.0I, N, N)
    m.p_c = zeros(N)
    m.p_σ = zeros(N)
    m.σ = σ_init
    h_σ!(m, gen)
    eigendecomp!(m)
    sqrtC!(m)
    invsqrtC!(m)
    m
  end
end


#-------------------------
#  Noise, including ShapedNoise and SphericalNoise

abstract type Noise end

mutable struct SphericalNoise <: Noise    # aka E
  values::Population
  weightedAvg::Vector
  flattened::Bool
  SphericalNoise(members) =
      (sn = new(); Noise(sn, members))
  SphericalNoise(members, structure::SortStructure) =
      (sn = new(); Noise(sn, members, structure))
  SphericalNoise(members, structure::SortStructure, model::CMAES_Model) =
      (sn = new(); Noise(sn, members, structure, model))
end

mutable struct ShapedNoise <: Noise     # aka W
  values::Population
  weightedAvg::Vector
  flattened::Bool
  ShapedNoise(members) =
      (sn = new(); Noise(sn, members))
  ShapedNoise(members, structure::SortStructure) =
      (sn = new(); Noise(sn, members, structure))
  ShapedNoise(members, structure::SortStructure, model::CMAES_Model) =
      (sn = new(); Noise(sn, members, structure, model))
end


#-------------------------
#  Noise, including ShapedNoise and SphericalNoise

struct Selection_Parms
  N::Integer
  μ::Integer
  λ::Integer
  direction::Symbol       # can be :min or :max
  loaded::Symbol
  beginBinding::Symbol
  segmentCount::Int
  includeCenter::Bool
  elitism::Bool
  η::Int64
  function Selection_Parms(N, μ, λ, direction, loaded, beginBinding, segmentCount, includeCenter, elitism, η)
    @assert μ < λ "Offspring population must be larger then parent population"
    new(N, μ, λ, direction, loaded, beginBinding, segmentCount, includeCenter, elitism, η)
  end
end

#-----------------------
# Selection Source Monitoring

abstract type SelectionSourceParms end

struct SingleSourceParms <: SelectionSourceParms
    sourceValues::Vector{Symbol}
    selectionType::Type
end

struct MultiSourceParms <: SelectionSourceParms
    sourceValues::Vector{Tuple{Symbol,Symbol}}
    selectionType::Type
end

abstract type SelectionSource end


mutable struct EliteSource <: SelectionSource
  source::Vector{Tuple{Symbol,Symbol}}
  fitness::Vector
  rank::Vector{Float64}
end

mutable struct CenterSource <: SelectionSource
  source::Vector{Tuple{Symbol,Symbol}}
  fitness::Vector
  rank::Vector{Float64}
end

mutable struct EliteCenterSource <: SelectionSource
  source::Vector{Tuple{Symbol,Symbol}}
  fitness::Vector
  rank::Vector{Float64}
end



const TripleSource = EliteCenterSource
DualSource = Union{CenterSource, EliteSource}
struct SingleSource <: SelectionSource
    source::Vector{Symbol}
    fitness::Vector
    rank::Vector{Float64}
end


#-------------------------
#  CMAES_System, including CMA_ES_System and IPop_System

struct CMAES_System <: System
  maxEvals::Int
  rParms::Reproduction_Parms  # constants used during  reproduction (aka Model_Parms)
  sParms::Selection_Parms     # constants used during selection
end


#-------------------------
#  CMAES_State, including CMA_ES_State and IPop_State

mutable struct CMAES_State <: State
   gen::Integer
   gen_shadow::Integer		    #shadow shadow
   evalCount::Integer
   evalCount_shadow::Integer		#shadow shadow
   sys::CMAES_System           # constants used during the running of the system (includes reproduction and selection parameters)
   pModel::CMAES_Model         # model from previous gen, used at the beginning of reproduction step
   pModel_shadow::CMAES_Model		# shadowged model from previous gen
   pOffspring::Population      # popnsize = μ    offspring from previous gen  (used for elitism)
   pOffspring_shadow::Population		# popnsize = μ	  shadowged offspring from previous gen (used for elitism)
   pE::SphericalNoise          # popnsize = μ    E from previous gen         (currenlty used for elitism) shadowGED doesnt matter nE is all same
   pW::ShapedNoise             # popnsize = μ    W from previous gen         (used for elitism)
   pW_shadow::ShapedNoise         # popnsize = μ    shadowged W from previous gen         (used for elitism)
   nE::SphericalNoise          # popnsize = λ    generated noise used to create samples
   nW::ShapedNoise             # popnsize = λ    shaped noise used to create samples
   nW_shadow::ShapedNoise         # popnsize = λ    shadowged shaped noise used to create samples
   nOffspring::Population      # popnsize = λ    samples used to create new model after selection
   nOffspring_shadow::Population  # popnsize = λ    shadowged samples used to create new model after selection
   sOffspring::Population      # popnsize = μ;   sorted and truncated (selection) changing sample distribution for model updating
   sOffspring_shadow::Population      # popnsize = μ;   shadowged sorted and truncated (selection) changing sample distribution for model updating
   sW::ShapedNoise             # popnsize = μ;   sorted and truncated (selection) changing sample distribution for model updating
   sW_shadow::ShapedNoise             # popnsize = μ;   shadowged sorted and truncated (selection) changing sample distribution for model updating
   nModel::CMAES_Model         # model after the reproductive step complete
   nModel_shadow::CMAES_Model     # shadowged model after the reproductive step complete
   status::Symbol
   status_shadow::Symbol			# shadowged status
   best::Tuple                  # shadowged best chromosome
   best_shadow::Tuple

  function CMAES_State(model::CMAES_Model, sys::CMAES_System, f::RealFitness, runInfo = NoMonitor(), verbose = false)
    state = new()
    setup!(state, model, sys, f)
    evolve!(state, f, runInfo, verbose)
    state
  end

  function CMAES_State(model::CMAES_Model, sys::CMAES_System, f::RealFitness, restart::RestartState,
                       runInfo = NoMonitor(), verbose = false)
    state = new()
    setup!(state, model, sys, f)
    evolve!(state, f, restart, runInfo, verbose)
    state
  end
end


#-------------------------
#  CMAES Restart

abstract type CMAES_Restart <: Restart end

mutable struct CMAES_RestartBase <: CMAES_Restart
  ignoreStagnation::Bool
  historyWindow::Union{Symbol, Int}
  η::Float64
  η_μ::Float64
  η_λ::Float64
  tol_f::Float64          # const   CMAES only (probably)
end

mutable struct CMAES_RestartFull <: CMAES_Restart
  ignoreStagnation::Bool
  historyWindow::Int
  η::Float64
  η_μ::Float64
  η_λ::Float64
  initλ::Int
  tol_f::Float64          # CMAES only (probably)
  tol_x::Float64          # CMAES only
  tol_x_::Float64          #uses σ from shadow model to calculate tol_x_
  center_init::Vector
  σ_init::Float64
end

#-------------------------
#  CMAES RunInfo

mutable struct CMAES_RunInfo <: RunInfo
  data::Dict
  uq::Tuple{Int64,Float64}
  med::Tuple{Int64,Float64}
  lq::Tuple{Int64,Float64}
  sourceValues::SelectionSourceParms
  allSourceValues::MultiSourceParms
  function CMAES_RunInfo(sys::CMAES_System)
    ri = new()
    ri.data = Dict()
    update!(ri, sys)
    ri
  end
end

abstract type DE_Parents end

struct DirectParents <: DE_Parents
  diffIndv1::Int64
  diffIndv2::Int64
  target::Int64
end

struct ClassicParents <: DE_Parents
  diffIndv1::Int64
  diffIndv2::Int64
  base::Int64
  target::Int64
end

struct BestParents <: DE_Parents
  diffIndv1::Int64
  diffIndv2::Int64
  best::Int64
  target::Int64
end

# often an individual if # of offpsring = 1, 
#    or a vector of indidividuals if there is more than one per target
struct DE_Individual
  genes::Vector         
end

mutable struct SlotSamplers
  samplers::Dict
  function SlotSamplers()
    s = new()
    s.samplers = Dict()
    s
  end
end

mutable struct PopnInitializer
  initializers::Dict
  function PopnInitializer()
    p = new()
    p.initializers = Dict()
    p
  end
end

struct DE_SamplingParms
  CR::Float64
  diffWeight::Float64
  donorSel::Symbol 
  donorSlope::Float64
  parentSel::Symbol 
  parentSlope::Float64
  localDonor::Float64
  template::Vector{Symbol}
end

struct DE_System <: System
  maxEvals::Int64
  popnSize::Int64
  chrLength::Int64    # for DE chromosome length = problem dimension
  selection::Function
  segCount::Int64
  randomizeSlots::Bool
  initialize::Function
  rParms::DE_SamplingParms        # sampling is actually reproduction in DE (selection and crossover) - so labled rParms
  samplers::Array{Function}
  evalsPerGen::Int64
end

const DE_Fitness = BoundedFitness{Float64}

function DE_Fitness(objFn::Function, dimension::Int64, direction::Symbol, 
                    optimalValue::Float64, epsilon::Float64, initMax::Float64, initMin::Float64)
  DE_Fitness(objFn, dimension, direction, optimalValue, epsilon, fill(initMax, dimension), fill(initMin, dimension))
end

objfn(f::Fitness) = f.objFn
returntype(f::DE_Fitness) = Float64

mutable struct Contenders
  targets::Population
  samples::Array{RegularPopulation}
  alternates::Array{RegularPopulation}
end

mutable struct DE_State <: State
  gen::Int64
  evalCount::Int64
  contenders::Contenders
  popn::Population
  status::Symbol
  best::Tuple
  sys::DE_System
  function DE_State(sys::DE_System, fit::DE_Fitness)
    state = new()
    setup!(state, sys, fit)
    state
  end
end

#-------------------------
#  DE Restart

abstract type DE_Restart <: Restart end

struct DE_RestartBase <: DE_Restart
  ignoreStagnation::Bool
  historyWindow::Union{Symbol, Int}
  η::Int
  tol_f::Float64          # const   CMAES only (probably)
end

struct DE_RestartFull <: DE_Restart
  ignoreStagnation::Bool
  historyWindow::Int
  η::Int
  tol_f::Float64
  initPopnSize::Int
end

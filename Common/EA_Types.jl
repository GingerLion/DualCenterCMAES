abstract type Monitor end
abstract type System end
abstract type State end
abstract type Restart end
abstract type Fitness end
abstract type Population end
abstract type ReturnInfo end
abstract type ReReturnInfo <: ReturnInfo end


const List = Union{DataStructures.Cons, DataStructures.Nil, Vector}

#-------------------------
#  EA_Fitness

struct UnBoundedFitness{T} <: Fitness
  objFn::Function
  dimension::Int64
  direction::Symbol     # either :min or :max respectively when minimizing or maximizing the function
  optimalValue::T       # if known, else set to Inf or -Inf depending on direction
  optimum::Vector       # if known, else set to [NaN, ..., NaN]
  epsilon::T            # difference from goal that is considered a success
end

struct BoundedFitness{T} <: Fitness
  objFn::Function
  dimension::Int64
  direction::Symbol     # either :min or :max respectively when minimizing or maximizing the function
  optimalValue::T       # if know else set to Inf or -Inf depending on direction
  optimum::Vector       # if known, else set to [NaN, ..., NaN]
  epsilon::T            # difference from goal that is considered a success
  initMax::Vector     # upper bound on the solution vector
  initMin::Vector     # lower bound on the solution vector
end

#-------------------------
#  EA_Population

mutable struct RegularPopulation <: Population
  members
  fitness
  direction::Symbol
  evaluated::Bool
end

mutable struct SortedPopulation <: Population
  members                   # holds the members from the original population - should not be modified
  fitness
  sortOrder::Vector       # sort for members [1:λ], used for reseting index in a truncation where μ_new > μ_old
  index::Vector           # sort for members [1:μ], used for indexing into the members
  lt::Function
  direction::Symbol
end

mutable struct SortStructure      # used by other populations to match sort order of the SortedPopulation
  fitness
  sortOrder::Vector       # sort for members [1:λ], used for reseting μ in a truncation where μ_new > μ_old
  index::Vector           # sort for members [1:μ], used for indexing into the members
                          #   -  naturally causes truncated out of bound errors
  lt::Function
  direction::Symbol
end

## created from the population, not for the population,
##    which is based on a 2d matrix of real values
struct Member
  chromosome
  fitness
end

#-------------------------
#  RunInfo and ReturnInfo

abstract type RunInfo <: Monitor end

mutable struct EA_RunInfo <: RunInfo
  data::Dict
  EA_RunInfo() = new(Dict())
end

struct NoMonitor <: Monitor end

const Storage = RunInfo           # depreciated
const NoStorage = NoMonitor       # depreciated

mutable struct EA_ReturnInfo  <: ReturnInfo
  state::State
  runInfo::Monitor
end

mutable struct EA_ReReturnInfo <: ReReturnInfo
  state::List
  runInfo::Monitor
  best::Tuple         # three tuple (chr, fit, rep) where rep is the repitition/restart that found the best
  best_shadow::Tuple
  best_overall::Tuple
  best_overall_::Tuple
  popnSize::List      # list of (mu, lambda) pairs for CMAES
  localEvals::List
  localEvals_::List
  totalEvals::List
  totalEvals_::List
  stagnation::List
  bad_count::List
  good_count::List
end


#-------------------------
#  RestartState (also includes BestFitHistory, which is used to check for a restart condition)

mutable struct BestFitHistory
  history::Vector{Float64}
  ptr::Int
  windowSize::Int
  direction::Symbol
  function BestFitHistory(windowSize::Int, direction::Symbol)
    bfh= new()
    bfh.windowSize = windowSize
    bfh.history = fill(Inf, windowSize)
  bfh.ptr = 1
  bfh.direction = direction
  bfh
  end
end

mutable struct BestChrHistory
  history::List
  ptr::Int
  windowSize::Int
  direction::Symbol
  function BestChrHistory(windowSize::Int, direction::Symbol)
    bch= new()
    bch.windowSize = windowSize
    bch.history = nil()
  bch.ptr = 1
  bch.direction = direction
  bch
  end
end

mutable struct RestartState
  rep::Int
  totalEvals::Integer
  totalEvals_::Integer
  bfHist::BestFitHistory
  bfHist_::BestFitHistory
  bcHist_::BestChrHistory
  bcHist_fitnesses_::BestFitHistory
  sortOrder::Vector{Int}
  stagnation::Vector{Bool}
  stagnation_::Vector{Bool}
  shouldRestart::Bool
  stagReason::Array{String}
  stagReason_::Array{String}
  parms::Restart

  function RestartState(bfHist::BestFitHistory, bcHist_::BestChrHistory, bcHist_fitnesses_::BestFitHistory, parms::Restart)
    restart = new()
    restart.rep = 0
    restart.totalEvals = 0
    restart.totalEvals_ = 0
    restart.bfHist = bfHist
    restart.bfHist_ = deepcopy(bfHist)
    restart.bcHist_ = bcHist_
    restart.bcHist_fitnesses_ = bcHist_fitnesses_
    restart.sortOrder = Vector{Int}(undef,length(bcHist_.windowSize))
    restart.shouldRestart = false
    restart.stagReason = Array{String, 1}()
    restart.stagReason_ = Array{String, 1}()
    restart.parms = parms
    restart
  end
end

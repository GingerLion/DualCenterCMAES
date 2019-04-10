#----------------------------------------------------
# mutable struct DE_State <: State
#----------------------------------------------------
#   gen::Int64
#   evalCount::Int64
#   contenders::Contenders
#   popn::Population
#   status::Symbol
#   best::Tuple
#   sys::DE_System
#----------------------------------------------------
#----------------------------------------------------
# Inner Constructor
#----------------------------------------------------
#    DE_State(sys::DE_System, fit::DE_Fitness)
#----------------------------------------------------

# needed for restarting since CMAES evolves the state and evolve! uses runInfo
function DE_State(sys::DE_System, fit::DE_Fitness, runInfo::Monitor, verbose::Verbose)
  state = DE_State(sys, fit)
  
  if atlevel(verbose, RunLevel)  
  	println(state) 
  end
  
  if monitorable(runInfo) 
  	monitor!(runInfo, state, restart, fit) 
  end
  state
end

function DE_State(sys::DE_System, fit::DE_Fitness, restart::RestartState, runInfo::Monitor, verbose::Verbose)
  DE_State(sys, fit, runInfo, verbose)
end

function setup!(state::DE_State, sys::DE_System, fit::DE_Fitness)
  state.sys = sys
  state.gen = 0
  evolvable!(state)
  state.popn = RegularPopulation(sys.initialize(sys), fit)
  state.evalCount = sys.evalsPerGen
  best!(state)
  found!(state, fit)        # best!() must be run before found!()
end

evolvable(state::DE_State, restart::RestartState) = evolvable(state)
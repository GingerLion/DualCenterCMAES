
#-----------------------------------------------------------------------
#  Restart System
#-----------------------------------------------------------------------
# struct CMAES_Restart <: Restart
#   testCount::Int
#   η::Int, 
#   η_μ::Int, 
#   η_λ::Int,
#   initλ::Int,
#   tol_f::Float64          # const   CMAES only (probably)
#   tol_x::Float64          # const   CMAES only
#-----------------------------------------------------------------------
# mutable struct RestartState
#   rep::Int
#   totalEvals::Integer
#   bfHist::BestFitHistory
#   stagnation::Vector{Bool}
#   parms::Restart
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Outer Constructors
#-----------------------------------------------------------------------

function DE_Restart(; historyWindow = :default, ignoreStagnation = false, η = 2, tol_f = 10.0^-12)
  DE_RestartBase(ignoreStagnation, historyWindow, η, tol_f)
end

function DE_RestartFull(rsys::DE_RestartBase, state::DE_State)
  initPopnSize = popnsize(state)
  hw = ((rsys.historyWindow == :default) ? 10 + ceil(Int, 30 * dimensions(state) / popnsize(state)) 
                                         : rsys.historyWindow)
  DE_RestartFull(rsys.ignoreStagnation, hw, rsys.η, rsys.tol_f, initPopnSize)
end


#-----------------------------------------------------------------------
# General Functions - RestartState only
#-----------------------------------------------------------------------

scalepopnsize(restart::RestartState) = restart.parms.η
initpopnsize(restart::RestartState) = restart.parms.initPopnSize
tol_f(restart::RestartState) = restart.parms.tol_f


#-----------------------------------------------------------------------
# RestartState & DE_State interactions
#-----------------------------------------------------------------------

function stagnationcriteria(state::DE_State, restart::RestartState)
#  (model, bfh, gen) = (currentmodel(state), bestfithist(restart), currentgen(state))
  [equalfunvalhist(state, restart)]
#   noeffectaxis(model, gen),
#   negeigerr(state),
#   complexeigerr(state),
#   zeroeigerr(state),
#   noeffectcoord(model)]
end

function nextpopnsize(sys::DE_System, restart::RestartState)
  popnsize(sys) * scalepopnsize(restart)
end

function next!(restart::RestartState, state::DE_State)
  restart.rep += 1
  histWindow = 10 * restart.rep + ceil(Int, 30 * dimensions(state) / initpopnsize(restart))
  restart.bfHist = BestFitHistory(histWindow, direction(state))
  restart.shouldRestart = false
end

#-----------------------------------------------------------------------
#  Restart Criteria
#-----------------------------------------------------------------------

# equalfunvalhist 
# The best fitness values of the last histWindow = 10 + ceil(30n/λ) generations 
# have a difference between their maximum and minimum values smaller than T_x
function equalfunvalhist(state::DE_State, restart::RestartState) 
  currentFit = copy(fitness(population(state)))
  equalfunvalhist(currentFit, restart)
end


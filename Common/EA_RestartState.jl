#--------------------------------------
# mutable struct RestartState
#--------------------------------------
#   rep::Int
#   totalEvals::Integer
#   bfHist::BestFitHistory
#   stagnation::Vector{Bool}
#	shouldRestart::Bool
#   parms::Restart
#--------------------------------------
# inferface Restart   (expected instance variables)
#   testCount::Int
#   ignoreStagnation::Bool
#   historyWindow::Union{Symbol, Int}
#   η::Int,
#--------------------------------------
#--------------------------------------
# Outer Constructor
#--------------------------------------

function RestartState(rsys::Restart, state::State)
  bfHist = BestFitHistory(rsys.historyWindow, direction(state))
  restart = RestartState(bfHist, rsys)
  stagnationupdate!(restart, state)
  restart
end

##--------------------------------------
##  General Functions

rep(restart::RestartState) = restart.rep
evals(restart::RestartState) = restart.totalEvals
evals_(restart::RestartState) = restart.totalEvals_
bestfithist(restart::RestartState) = restart.bfHist
bestfithist_(restart::RestartState) = restart.bfHist_
ignorestagnation(restart::RestartState) = restart.parms.ignoreStagnation
stagnation(restart::RestartState) = restart.stagnation
stagnflags(restart::RestartState) = restart.stagnation
shouldrestart(restart::RestartState) = restart.shouldRestart


function stagnationupdate!(restart::RestartState, state::State)
  if (currentgen(state) > 0)
  	restart.stagnation = stagnationcriteria(state, restart)
    restart.stagnation_ = stagnationcriteria_(state, restart)

    if any(restart.stagnation) && any(restart.stagnation_)
      restart.shouldRestart = true
    elseif any(restart.stagnation) && !any(restart.stagnation_)
      state.firstRequest = :normal
    elseif any(restart.stagnation_) && !any(restart.stagnation)
      state.firstRequest = :dualcenter
    end

  end
end

function update!(restart::RestartState, state::State, sys::System)
  restart.totalEvals += evalsPerGen(sys)
  if status_(state) != :found
      restart.totalEvals_ += evalsPerGen_(sys)
  else
      restart.totalEvals_ = restart.totalEvals_
  end
  restart.bfHist[] = bestfitness(best(state))
  restart.bfHist_[] = bestfitness(best_(state))
end

##--------------------------------------
##  General Restart Criteria

function equalfunvalhist(currentFit::Vector, restart::RestartState)
	historyFit = history(bestfithist(restart))
	allFit = vcat(currentFit, historyFit)
	zero_hist = (maximum(historyFit) - minimum(historyFit) == 0.0)
	tol_all = (maximum(allFit) - minimum(allFit) < tol_f(restart))
	zero_hist || tol_all
end

function equalfunvalhist_(currentFit::Vector, restart::RestartState)
	historyFit = history(bestfithist_(restart))
	allFit = vcat(currentFit, historyFit)
	zero_hist = (maximum(historyFit) - minimum(historyFit) == 0.0)
	tol_all = (maximum(allFit) - minimum(allFit) < tol_f(restart))
	zero_hist || tol_all
end

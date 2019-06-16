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
  bcHist_ = BestChrHistory(10, direction(state))
  restart = RestartState(bfHist, bcHist_, rsys)
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
  if (currentgen(state) > 0) || (currentgen_(state) > 0)
  	restart.stagnation = stagnationcriteria(state, restart)
    restart.stagnation_ = stagnationcriteria_(state, restart)

    #=if equalfunvalhist(state, restart)
      state.firstRequest = :normal
      restart.shouldRestart = true
    elseif equalfunvalhist_(state, restart)
      state.firstRequest = :dualcenter
      restart.shouldRestart = true=#
    if any(restart.stagnation) && any(restart.stagnation_)
        restart.shouldRestart = true
    elseif found(state) && (any(restart.stagnation_) || status_(state) == :stop)
        restart.shouldRestart = true
    elseif found_(state) && (any(restart.stagnation) || status(state) == :stop)
        restart.shouldRestart = true
    elseif status(state) == :stop && status_(state) == :stop
        restart.shouldRestart = true
    elseif any(restart.stagnation) && !any(restart.stagnation_)
      state.status = :stop
      state.firstRequest = :normal
      state.stopEvals = true
    elseif any(restart.stagnation_) && !any(restart.stagnation)
      state.status_shadow = :stop
      state.firstRequest = :dualcenter
      state.stopEvals_ = true
    end

  end
end

function update!(restart::RestartState, state::State, sys::System)
  if !stopEvals(state) restart.totalEvals += evalsPerGen(sys) end
  if !stopEvals_(state) restart.totalEvals_ += evalsPerGen_(sys) end
  restart.bfHist[] = bestfitness(best(state))
  restart.bfHist_[] = bestfitness(best_(state))
  restart.bcHist_[] = bestchromosome(best_(state)) #note that this is a list with its own setindex! function
  #println("bcHist_ = $(history(restart.bcHist_))")
  len = length(history(restart.bcHist_))

  w = Weights(normalize(map((i)->(log(len+chrlength(state)) - log(i)), 1:len), 1))
  center!_(state,  sum(map((x) -> w[x] * history(restart.bcHist_)[x], 1:length(w))))
  #center!_(state, mean(history(restart.bcHist_), w, dims=1))
  if len >=  length(restart.bcHist_) restart.bcHist_.history = nil() end
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

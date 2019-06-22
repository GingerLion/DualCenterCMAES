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
    if (found(state) || status(state) == :max_evals) restart.stagnation = fill(false, length(restart.stagnation)) end
    if (found_(state) || status_(state) == :max_evals) restart.stagnation_ = fill(false, length(restart.stagnation_)) end

    if any(restart.stagnation) && any(restart.stagnation_)
        restart.shouldRestart = true
        #println("case 1 shouldRestart")
    elseif found(state) && (any(restart.stagnation_) || status_(state) == :stop)
        restart.shouldRestart = true
        #println("case 2 shouldRestart")
    elseif found_(state) && (any(restart.stagnation) || status(state) == :stop)
        restart.shouldRestart = true
        #println("case 3 shouldRestart")
    elseif status(state) == :stop && status_(state) == :stop
        restart.shouldRestart = true
        #println("case 4 shouldRestart")
    elseif status(state) == :stop && status_(state) == :max_evals
        restart.shouldRestart = true
        #println("case 5 shouldRestart")
    elseif status_(state) == :stop && status(state) == :max_evals
        restart.shouldRestart = true
        #println("case 6 shouldRestart")
    elseif any(restart.stagnation) && !any(restart.stagnation_)
      state.status = :stop
      if status_(state) == :stop restart.shouldRestart = true end
      state.firstRequest = :normal
      state.stopEvals = true
      #println("case 7")
    elseif any(restart.stagnation_) && !any(restart.stagnation)
      state.status_shadow = :stop
      if status(state) == :stop restart.shouldRestart = true end
      state.firstRequest = :dualcenter
      state.stopEvals_ = true
      #println("case 8")
    end
  end
end

function update!(restart::RestartState, state::State, sys::System)
  if !stopEvals(state) restart.totalEvals += evalsPerGen(sys) end
  if !stopEvals_(state) restart.totalEvals_ += evalsPerGen_(sys) end

  if evolvable(state) restart.bfHist[] = bestfitness(best(state)) end
  if evolvable_(state)
      restart.bfHist_[] = bestfitness(best_(state))
      restart.bcHist_[] = bestchromosome(best_(state)) #note that this is a list with its own setindex! function
  end
  #println("bcHist_ = $(history(restart.bcHist_))")
  len = length(history(restart.bcHist_))
  if len > 0
      w = Weights(normalize(map((i)->(log(len+chrlength(state)) - log(i)), 1:len), 1))
      center!_(state,  sum(map((x) -> w[x] * history(restart.bcHist_)[x], 1:length(w))))
      #center!_(state, mean(history(restart.bcHist_), w, dims=1))
  else
      try
          center!_(state, bestfitness(best_(state)))
      catch
          center!_(state, center_(state))
      end
  end
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

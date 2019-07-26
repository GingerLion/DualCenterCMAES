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
  bcHist_fitnesses_ = BestFitHistory(10, direction(state))
  restart = RestartState(bfHist, bcHist_, bcHist_fitnesses_, rsys)
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

function stagReason!(restart::RestartState, msg::Array{String})
    for i=1:length(restart.stagnation)
        if restart.stagnation[i]
            if !any(j -> (j == msg[i]), restart.stagReason)
                restart.stagReason = vcat(restart.stagReason, msg[i])
            end
        end
    end
end

function stagReason!_(restart::RestartState, msg::Array{String})
    for i=1:length(restart.stagnation_)
        if restart.stagnation_[i]
            if !any(j -> (j == msg[i]), restart.stagReason_)
                restart.stagReason_ = vcat(restart.stagReason_, msg[i])
            end
        end
    end
end

function stagnationupdate!(restart::RestartState, state::State)
  if (currentgen(state) > 0) || (currentgen_(state) > 0)
  	restart.stagnation = stagnationcriteria(state, restart)
    restart.stagnation_ = stagnationcriteria_(state, restart)
    if (found(state) || status(state) == :max_evals) restart.stagnation = fill(false, length(restart.stagnation)) end
    if (found_(state) || status_(state) == :max_evals) restart.stagnation_ = fill(false, length(restart.stagnation_)) end

    msg = ["hist", "tolx", "zeroAxis", "negEigVal", "complexEigVal", "zeroEigVal", "zeroCoord", "cond"]

    if any(restart.stagnation) && any(restart.stagnation_)
        stagReason!(restart, msg)
        stagReason!_(restart, msg)
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
      stagReason!(restart, msg)
      state.status = :stop
      if status_(state) == :stop restart.shouldRestart = true end
      state.firstRequest = :normal
      state.stopEvals = true
      #println("case 7")
    elseif any(restart.stagnation_) && !any(restart.stagnation)
      stagReason!_(restart, msg)
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
      if !isWindowFull(restart)
          restart.bcHist_[] = bestchromosome(best_(state)) #note that this is a list with its own setindex! function
          #println("Chromosome Window => $(history(restart.bcHist_))")
          restart.bcHist_fitnesses_[] = bestfitness(best_(state))
          #println("Fitness Window => $(history(restart.bcHist_fitnesses_))")
          #println("Window is not full. Length of window is $(length(history(restart.bcHist_))).")
      elseif isWindowFull(restart) && canEnterWindow(deepcopy(bestfitness(best_(state))), history(restart.bcHist_fitnesses_), direction(restart.bcHist_))
          #println("-----------------------------------------\n Can enter window.")
          tmp = deepcopy(history(restart.bcHist_))
          #println("tmp before delete => $(tmp)")
          deleteat!(tmp,length(restart.bcHist_)) #delete at 10 (oldest member of list)
          #println("tmp after delete => $(tmp)")
          #println("fitness window before delete => $(history(restart.bcHist_fitnesses_))")
          deleteat!(restart.bcHist_fitnesses_.history, 1) # delete at 1 (oldest fitness of the vector)
          #println("fitness window before delete => $(history(restart.bcHist_fitnesses_))")
          restart.bcHist_.history = nil() #empty list
          #dont have to empty fitness vector
          for i in Iterators.reverse(tmp) #put old values back from tmp to the List
            restart.bcHist_.history = cons(i, restart.bcHist_.history) # List
          end
          restart.bcHist_.history = cons(bestchromosome(best_(state)), restart.bcHist_.history) # add new chromosome (most recent and head of list)
          #println("Chromosome Window with new chromosome => $(history(restart.bcHist_))")
          push!(restart.bcHist_fitnesses_.history, bestfitness(best_(state))) # append new fitness of new chromosome (most recent at end of vector)
          #println("fitness window after new fitness => $(history(restart.bcHist_fitnesses_))\n-----------------------------------------")
      else
          println("Couldn't enter window.")
      end
  end
  # NON-ELITIST window code
  #=restart.bcHist_[] = bestchromosome(best_(state))
  len = length(history(restart.bcHist_))
  if len >  length(restart.bcHist_)
    tmp = deepcopy(history(restart.bcHist_)) # List
    deleteat!(tmp,len) # List
    restart.bcHist_.history = nil()
    for i in Iterators.reverse(tmp)
      restart.bcHist_.history = cons(i, restart.bcHist_.history) # List
    end
end=# #NON - ELITIST Window code
  len = length(history(restart.bcHist_))
  if len > 1
      w = Weights(normalize(map((i)->(log(len + 15 + (4 * log(chrlength(state)))) - log(i)), 1:len), 1))
      #println("lenw = $(length(w)) len = $(len), \n w = $(w), history = $(history(restart.bcHist_))")
      center!_(state,  sum(map((x) -> w[x] * history(restart.bcHist_)[x], 1:length(w))))
      #center!_(state, mean(history(restart.bcHist_), w, dims=1))
  else
      try
          center!_(state, bestchromosome(best_(state)))
      catch
          center!_(state, center_(state))
      end
  end

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

isWindowFull(restart::RestartState) = (length(history(restart.bcHist_)) == length(restart.bcHist_)) ? true : false

function canEnterWindow(newFit::Float64, historyWindow::Vector{Float64}, direction::Symbol)
    canEnter = false
    if direction == :min
        for i in historyWindow
            if newFit < i
                canEnter = true
            end
        end
    elseif direction == :max
        for i in historyWindow
            if newFit > i
                canEnter = true
            end
        end
    end
    return canEnter
end

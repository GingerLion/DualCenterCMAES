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

function update!(restart::RestartState, state::State, sys::System, f::Fitness)
  if !stopEvals(state) restart.totalEvals += evalsPerGen(sys) end
  if !stopEvals_(state) restart.totalEvals_ += evalsPerGen_(sys) end

  if evolvable(state) restart.bfHist[] = bestfitness(best(state)) end
  if evolvable_(state)
      restart.bfHist_[] = bestfitness(best_(state))
      #updateSlidingWindow(state, restart, f)
      updateEliteWindowUnSorted(state, restart, f)
      #updateEliteWindowSorted(state, restart, f)
  end
  #println("normal evals = ",evals(restart),", dc evals = ",evals_(restart))
end

##--------------------------------------
##  General Restart Criteria

function equalfunvalhist(currentFit::Vector, restart::RestartState)
	historyFit = history(bestfithist(restart))
	allFit = vcat(currentFit, historyFit)
	zero_hist = (maximum(historyFit) - minimum(historyFit) == 0.0)
	tol_all = (maximum(allFit) - minimum(allFit) < tol_f(restart))
    #if (zero_hist || tol_all) println("<HIST> triggered") end
	zero_hist || tol_all
end

function equalfunvalhist_(currentFit::Vector, restart::RestartState)
	historyFit = history(bestfithist_(restart))
	allFit = vcat(currentFit, historyFit)
	zero_hist = (maximum(historyFit) - minimum(historyFit) == 0.0)
	tol_all = (maximum(allFit) - minimum(allFit) < tol_f(restart))
    #if (zero_hist || tol_all) println("<HIST> triggered") end
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
# constantly moving sliding window, solutions will only last windowSize generations
function updateSlidingWindow(state::State, restart::RestartState, f::Fitness)
    restart.bcHist_[] = bestchromosome(best_(state))
    len = length(history(restart.bcHist_))
    if len >  length(restart.bcHist_)
      tmp = deepcopy(history(restart.bcHist_)) # List
      #println("tmp chr window before => $(tmp)")
      deleteat!(tmp,len) # List
      #println("tmp chr window after delete=> $(tmp)")
      restart.bcHist_.history = nil()
      for i in Iterators.reverse(tmp)
          restart.bcHist_.history = cons(i, restart.bcHist_.history) # List
      end
    end
    #println("real chr window before update=> $(history(restart.bcHist_))")
    # update center_
    updateCenter!_(state, restart, history(restart.bcHist_), f)
end

#keeps the best solutions at all times but doenst sort them
function updateEliteWindowUnSorted(state::State, restart::RestartState, f::Fitness)
    if !isWindowFull(restart)
        restart.bcHist_[] = bestchromosome(best_(state)) #note that this is a list with its own setindex! function
        #println("Chromosome Window => $(history(restart.bcHist_))\n")
        restart.bcHist_fitnesses_[] = bestfitness(best_(state))
        #println("Fitness Window => $(history(restart.bcHist_fitnesses_))\n")
        #println("Window is not full. Length of window is $(length(history(restart.bcHist_)))\n")
    elseif isWindowFull(restart) && canEnterWindow(deepcopy(bestfitness(best_(state))), history(restart.bcHist_fitnesses_), direction(restart.bcHist_))
        #println("-----------------------------------------\n Can enter window.")
        tmp = deepcopy(history(restart.bcHist_))
        #println("tmp before delete => $(tmp)\n")
        deleteat!(tmp,length(restart.bcHist_)) #delete at 10 (oldest member of list)
        #println("tmp after delete => $(tmp)\n")
        #println("fitness window before delete => $(history(restart.bcHist_fitnesses_))\n")
        deleteat!(restart.bcHist_fitnesses_.history, 1) # delete at 1 (oldest fitness of the vector)
        #println("fitness window after delete => $(history(restart.bcHist_fitnesses_))\n")
        restart.bcHist_.history = nil() #empty list
        #dont have to empty fitness vector
        for i in Iterators.reverse(tmp) #put old values back from tmp to the List
          restart.bcHist_.history = cons(i, restart.bcHist_.history) # List
        end
        restart.bcHist_.history = cons(bestchromosome(best_(state)), restart.bcHist_.history) # add new chromosome (most recent and head of list)
        #println("Chromosome Window with new chromosome => $(history(restart.bcHist_))\n")
        push!(restart.bcHist_fitnesses_.history, bestfitness(best_(state))) # append new fitness of new chromosome (most recent at end of vector)
        #println("fitness window after new fitness => $(history(restart.bcHist_fitnesses_))\n-----------------------------------------")
        #println("Couldn't enter window.")
    end
    #update center_
    updateCenter!_(state, restart, history(restart.bcHist_), f)
end
#keeps the best solutions and also keeps the sorted
function updateEliteWindowSorted(state::State, restart::RestartState, f::Fitness)
    if !isWindowFull(restart)
        #println("Window not full.")
        restart.bcHist_[] = bestchromosome(best_(state))
        #println("Chromosome Window => $(history(restart.bcHist_))")
        restart.bcHist_fitnesses_[] = bestfitness(best_(state))
        #println("Fitness Window => $(history(restart.bcHist_fitnesses_))")
        restart.sortOrder = sortperm(history(restart.bcHist_fitnesses_))
        #println("Sort Order => $(restart.sortOrder)")
    elseif isWindowFull(restart) && canEnterWindow(deepcopy(bestfitness(best_(state))), history(restart.bcHist_fitnesses_), direction(restart.bcHist_))
        #println("--------------------------------------------------")
        # delete weakest fitness value
        #println("Sort Order => $(restart.sortOrder)")
        #println("Fitness Window Before Delete=> $(history(restart.bcHist_fitnesses_))")
        deleteat!(restart.bcHist_fitnesses_.history, restart.sortOrder[end])
        #println("Fitness Window After Delete=> $(history(restart.bcHist_fitnesses_))")
        # push new fitness value to the end of array
        push!(restart.bcHist_fitnesses_.history, bestfitness(best_(state)))
        #println("Fitness Window After Adding New Fitness => $(history(restart.bcHist_fitnesses_))")
        # reverse and collect in temporary variable because list members are added the head of the list
        tmp = collect(reverse(history(restart.bcHist_)))
        #println("Real Chromosome Window = $(history(restart.bcHist_))")
        #println("Reversed Temporary Chromosome Window => $(tmp)")
        # delete chromosome associaed with weakest fitness value from the temporary array
        deleteat!(tmp, restart.sortOrder[end])
        #println("Reversed Temporary Chromosome Window After Delete => $(tmp)")
        # push new chromosome to the end of tmp
        push!(tmp, bestchromosome(best_(state)))
        #println("Reversed Temporary Chromosome Window After Adding New Chr => $(tmp)")
        #empty the list before adding back the chromosome list elements
        restart.bcHist_.history = nil()
        # add back the chromosomes from tmp to the chromosome window list (since tmp is already the reversed list i can simply cons() each element in order)
        for i in tmp
            restart.bcHist_.history = cons(i, restart.bcHist_.history)
        end
        #println("New Real Chromosome Window => $(history(restart.bcHist_))")
        # update the sortOrder since a new fitness values were added
        restart.sortOrder = sortperm(history(restart.bcHist_fitnesses_))
        #println("New Sort Order => $(restart.sortOrder)")

    end
    #prepare tmp center_ update
    len = length(history(restart.bcHist_))
    tmp = Array{Array{Float64,1},1}(undef, len)
    for i=1:len
        tmp[i] = collect(reverse(restart.bcHist_.history))[restart.sortOrder[i]]
    end
    #println("Chromosomes in order of fitness and ready to be flattened => $(tmp)")
    #println("--------------------------------------------------")
    #update center_
    updateCenter!_(state, restart, tmp, f)
end

function updateCenter!_(state::State, restart::RestartState, window::Array, f::Fitness)
    len = length(history(restart.bcHist_))
    if len > 1
        w = Weights(LinearAlgebra.normalize(map((i)->(log(15 + (4 * log(chrlength(state)))) - log(i)), 1:len), 1))
        center!_(state, sum(map((x)->w[x] * window[x],1:length(w))))
    else
        try
            center!_(state, bestchromosome(best_(state)))
        catch
            center!_(state, center_(state))
        end
    end
    #change_scales!(state, restart, f)
end

function change_scales!(state::State, restart::RestartState, f::Fitness)
    if minimizing(population_shadow(state))
        # if the window center (best center) is more fit than the main center then generate solutions from that center,
        # otherwise generate no solutions from the best center
        if fitness(centerpopn_(currentmodel_(state),f)) < fitness(centerpopn(currentmodel_(state),f))
            state.evalCount_shadow += 2
            restart.totalEvals_ += 2
            orig_scale!(currentmodel_(state),0.0)
            best_scale!(currentmodel_(state),2.0)
        else
            state.evalCount_shadow += 2
            restart.totalEvals_ += 2
            orig_scale!(currentmodel_(state),1.0)
            best_scale!(currentmodel_(state),1.0)
        end
        #=if (!(ranksum_best(state) == 0) && (ranksum_best(state) < ranksum_orig(state)))
            orig_scale!(currentmodel_(state), 1.0)
            best_scale!(currentmodel_(state), 1.0)
        else # >=
            orig_scale!(currentmodel_(state), 1.5)
            best_scale!(currentmodel_(state), 0.5)
        end=#
    else
        if fitness(centerpopn_(currentmodel_(state),f)) > fitness(centerpopn(currentmodel_(state),f))
            state.evalCount_shadow += 2
            restart.totalEvals_ += 2 # most important
            orig_scale!(currentmodel_(state),0.0)
            best_scale!(currentmodel_(state),2.0)
        else
            state.evalCount_shadow += 2
            restart.totalEvals_ += 2 # most important
            orig_scale!(currentmodel_(state),1.0)
            best_scale!(currentmodel_(state),1.0)
        end
        #=if (!(ranksum_best(state)==0) && (ranksum_best(state) > ranksum_orig(state)))
            orig_scale!(currentmodel_(state), 1.0)
            best_scale!(currentmodel_(state), 1.0)
        else # >=
            orig_scale!(currentmodel_(state), 1.5)
            best_scale!(currentmodel_(state), 0.5)
        end=#
    end
end

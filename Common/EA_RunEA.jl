function runEA(state::State, f::Fitness, runInfo::Monitor, returnInfo::ReturnInfo, verbose::Verbose)
  while evolvable(state)
    evolve!(state, f, runInfo, verbose)
  end

  decgen!(state)
  returning(returnInfo, verbose)
end

function restarting!(state::State, restart::RestartState, f::Fitness, runInfo::Monitor, verbose::Verbose)
  if atlevel(verbose, RestartLevel)
    println(state, restart)
  end

  next!(restart, state)
  sys = system(state)
  sys = typeof(sys)(sys, f, restart, verbose)
  update!(runInfo, sys)

  if hastocatchup(state) == :dualcenter || hastocatchup(state) == :normal || isMaxAndStop(state)
      typeof(state)(sys, f, restart, runInfo, verbose, new = false, cur_state = deepcopy(state))
  else
      typeof(state)(sys, f, restart, runInfo, verbose, new = true, cur_state = deepcopy(state))
  end
end

function runEA(state::State, restart::RestartState, f::Fitness,
               runInfo::Monitor, returnInfo::ReReturnInfo, verbose::Verbose)
  while (evolvable(state, restart) || hastocatchup(state) == :normal || evolvable_(state, restart) || hastocatchup(state) == :dualcenter || isMaxAndStop(state)) || (status(state) == :stop && status_(state) == :stop)
    if shouldrestart(restart)
      #println("restarting...")
      #if (status(state) == :stop && status_(state) == :stop) println("both systems = :stop") end
      update!(returnInfo, state, restart)
      state = restarting!(state, restart, f, runInfo, verbose)
    end

    #this block of code was commented out because now we're running both systems to max_evals a.k.a fixed budget
    #=if found(state) && !found_(state)
        system(state).maxEvals = first(runInfo[:total_evals])
    elseif found_(state) && !found(state)
        system(state).maxEvals = first(runInfo[:total_evals_])
    end=#

    #println("runEA: status = $(status(state)), status_shadow = $(status_(state)), maxEvals = $(system(state).maxEvals)")
    evolve!(state, f, restart, runInfo, verbose)
  end
  #println("hastocatchup = $(hastocatchup(state))")
  #println("bestfitOverall = $(bestfitoverall(state)) \n bestfitOverall_shadow = $(bestfitoverall_(state))")
  decgen!(state)
  decgen!_(state)
  update!(returnInfo, state, restart)
  returning(returnInfo, restart, verbose) # WHY?
end

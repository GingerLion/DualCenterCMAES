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
  typeof(state)(sys, f, restart, runInfo, verbose)
end

function runEA(state::State, restart::RestartState, f::Fitness,
               runInfo::Monitor, returnInfo::ReReturnInfo, verbose::Verbose)
  while (evolvable(state, restart) || evolvable_(state, restart))
    if shouldrestart(restart)
      update!(returnInfo, state, restart)
      state = restarting!(state, restart, f, runInfo, verbose)
    end

    evolve!(state, f, restart, runInfo, verbose)
  end

  decgen!(state)
  decgen!_(state)
  update!(returnInfo, state, restart)
  returning(returnInfo, restart, verbose) # WHY?
end

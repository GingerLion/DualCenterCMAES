#-----------------------------------------------------------------------
# evolve! functions
#-----------------------------------------------------------------------
# Error checks the evolvability of the population
# Calls evolvepopulation!
# Monitors the resulting state of the system
#-----------------------------------------------------------------------

function evolve!(state::State, f::Fitness; runInfo = NoMonitor(), verbose = NotVerbose())
  evolve!(state, f, runInfo, verbose)
end

function evolve!(state::State, f::Fitness, runInfo::Monitor, verbose::Verbose)
  if evolvable(state) || evolvable_(state)
    if currentgen_(state) == 0
    	   runInfo.sourceValues = SelectionSourceParms(system(state), currentmodel_(state))
    end
    if evolvable(state)
        evolvepopn!(state, f)
        gtmaxevals!(state)
    end
    if evolvable_(state)
        evolvepopn!_(state, f)
        gtmaxevals!_(state)
    end
    if monitorable(runInfo)        monitor!(runInfo, state, f) end
    if atlevel(verbose, RunLevel)  println(state) end
  end
  status(state)
end

#-----------------------------------------------------------------------
# evolve! functions with restart/population-resizing capability
#-----------------------------------------------------------------------

function evolve!(state::State, f::Fitness, restart::RestartState;
                 runInfo = NoMonitor(), verbose = NotVerbose())
  evolve!(state, f, restart, runInfo, verbose)
end

function evolve!(state::State, f::Fitness, restart::RestartState, runInfo::Monitor, verbose::Verbose)
  if (evolvable(state, restart) || status(state) == :stop) || (evolvable_(state, restart) || status_(state) == :stop)
    if currentgen_(state) == 0
      	runInfo.sourceValues = SelectionSourceParms(system(state), currentmodel_(state))
    end
    if evolvable(state, restart)
        evolvepopn!(state, f)
    end
    if evolvable_(state, restart)
        evolvepopn!_(state, f)
    end

    update!(restart, state, system(state))
    gtmaxevals!(state, restart)

    if (!ignorestagnation(restart) && ((evolvable(state, restart) || status(state) == :stop) || (evolvable_(state, restart) || status_(state) == :stop)))
      stagnationupdate!(restart, state)
    end


    if monitorable(runInfo)        monitor!(runInfo, state, restart, f) end
    if atlevel(verbose, RunLevel)  println(state) end
  end
  status(state)
end

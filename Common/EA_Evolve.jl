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
  if evolvable(state)
    evolvepopn!(state, f)

    if evolvable(state)
      gtmaxevals!(state)
    end

    if monitorable(runInfo)        monitor!(runInfo, state, f) end
    if atlevel(verbose, RunLevel)  println(state) end
  end
  if status(state) == :found
      status(state)
  elseif status_(state) == :found
      status_(state)
end

#-----------------------------------------------------------------------
# evolve! functions with restart/population-resizing capability
#-----------------------------------------------------------------------

function evolve!(state::State, f::Fitness, restart::RestartState;
                 runInfo = NoMonitor(), verbose = NotVerbose())
  evolve!(state, f, restart, runInfo, verbose)
end

function evolve!(state::State, f::Fitness, restart::RestartState, runInfo::Monitor, verbose::Verbose)
  if evolvable(state, restart)
    evolvepopn!(state, f)
    update!(restart, state, system(state))

    if evolvable(state, restart)
      gtmaxevals!(state, restart)
    end

    if (!ignorestagnation(restart) && evolvable(state, restart))
      stagnationupdate!(restart, state)
    end


    if monitorable(runInfo)        monitor!(runInfo, state, restart, f) end
    if atlevel(verbose, RunLevel)  println(state) end
  end

  status(state)
end

abstract type Verbose end
abstract type VerboseRun <: Verbose end
abstract type VerboseRestart <: VerboseRun end
abstract type VerboseReturn <: VerboseRestart end

struct ReturnLevel <: VerboseReturn end
struct RestartLevel <: VerboseRestart end
struct RunLevel <: VerboseRun end
struct NotVerbose <: Verbose end

isverbose(v::Verbose) = true
isverbose(nv::NotVerbose) = false

atlevel(verbose::Verbose, level::Type) = (supertype(level) <: supertype(typeof(verbose)))

#-----------------------------------------------
#  default EA println() routines used for verbose mode


function println(returnInfo::ReturnInfo)
	state = returnInfo.state
	(bestChr, bestFit) = best(state)
	(bestChr_bug, bestFit_bug) = best_(state)
	println("\n$(state.status): gen = $(gen(state)),  evals = $(evals(state))")
	println("best: fit[$bestFit]")
	println("best: bug_fit[$bestFit_bug]")
	println("best chr$bestChr")
	println("best bug chr$bestChr_bug")
end

function println(returnInfo::ReReturnInfo, restart::RestartState)
    state = finalstate(returnInfo)
    (bestChr, bestFit, bestRep) = returnInfo.best
	(bestChr_bug, bestFit_bug) = returnInfo.best_bug
    currentRep = reps(returnInfo)
    println("\n$(state.status): after $(evals(restart)) evals, at gen = $(gen(state)) within rep = $(reps(returnInfo))")
    print("best: ")
    if currentRep == bestRep
    	println("fit[$bestFit]")
		println("bug fit[$bestFit_bug]")
    else
    	println("found during rep = $bestRep")
    	println("      fit[$bestFit]")
		println("  bug fit[$bestFit_bug]")
    end
    println("      chr$bestChr")
	println("  bug chr$bestChr_bug")
end

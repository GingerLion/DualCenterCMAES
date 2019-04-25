#----------------------------------------------------------------------
# EA_ReturnInfo - assumed instance variables
#----------------------------------------------------------------------
#   state
#   runInfo
#----------------------------------------------------------------------

#--------------------------------------
# ReturnInfo - Constructors

ReturnInfo(state::State, runInfo::Monitor) = EA_ReturnInfo(state, runInfo)


#--------------------------------------
# ReReturnInfo - Updates

# None needed because State and Monitor will change inside ReturnInfo as it changes outside

#--------------------------------------
# EA_ReturnInfo - Control Function


function returning(returnInfo::ReturnInfo, verbose::Verbose)
  if atlevel(verbose, ReturnLevel)
  	println(returnInfo)
  end
  if monitorable(returnInfo.runInfo)
  	reverse!(returnInfo.runInfo)
  	collect!(returnInfo.runInfo)
  end
  returnInfo
end

runinfo(ri::ReturnInfo) = ri.runInfo
returnstate(ri::ReturnInfo) = ri.state

#----------------------------------------------------------------------
# mutable struct ReReturnInfo <: EA_ReturnInfo
#----------------------------------------------------------------------
#   state
#   runInfo
#   best
#   popnSize
#   localEvals
#   totalEvals
#   stagnation
#----------------------------------------------------------------------

#--------------------------------------
# ReReturnInfo - Constructors

function ReReturnInfo(state::State, runInfo::Monitor)
	bestMember = best(state)
    bestMember_shadow = best_(state)
    EA_ReReturnInfo(nil(), runInfo, (bestMember[1], bestMember[2], 1), (bestMember_shadow[1], bestMember_shadow[2], 1), nil(), nil(), nil(), nil(), nil(), nil())
end

#--------------------------------------
# ReReturnInfo - Updates

function update!(rinfo::ReReturnInfo, state::State, restart::RestartState)
  rinfo.state       = cons(deepcopy(state), rinfo.state)
  rinfo.popnSize    = cons(popnsize(state), rinfo.popnSize)
  rinfo.localEvals  = cons(evals(state), rinfo.localEvals)
  rinfo.localEvals_ = cons(evals_(state), rinfo.localEvals_)
  rinfo.totalEvals  = cons(evals(restart), rinfo.totalEvals)
  rinfo.totalEvals_  = cons(evals_(restart), rinfo.totalEvals_)
  rinfo.stagnation  = cons(stagnflags(restart), rinfo.stagnation)
  #restart.parms.σ_init *= 2
  if better(state, rinfo.best)
    rinfo.best = (bestchromosome(state), bestfitness(state), rep(restart))
  end
  if better_(state, rinfo.best_shadow)
    rinfo.best_shadow = (bestchromosome_(state), bestfitness_(state), rep(restart))
  end
end

function reverse!(rinfo::ReReturnInfo)
  rinfo.state = reverse(rinfo.state)
  rinfo.popnSize = reverse(rinfo.popnSize)
  rinfo.localEvals = reverse(rinfo.localEvals)
  rinfo.localEvals_ = reverse(rinfo.localEvals_)
  rinfo.totalEvals = reverse(rinfo.totalEvals)
  rinfo.totalEvals_  = reverse(rinfo.totalEvals_)
  rinfo.stagnation = reverse(rinfo.stagnation)
end

function collect!(rinfo::ReReturnInfo)
  rinfo.state = collect(rinfo.state)
  rinfo.popnSize = collect(rinfo.popnSize)
  rinfo.localEvals = collect(rinfo.localEvals)
  rinfo.localEvals_ = collect(rinfo.localEvals_)
  rinfo.totalEvals = collect(rinfo.totalEvals)
  rinfo.totalEvals_  = collect(rinfo.totalEvals_)
  rinfo.stagnation = collect(rinfo.stagnation)
end


#--------------------------------------
# ReReturnInfo - Control Functions

reps(returnInfo::ReReturnInfo) = length(returnInfo.state)
finalstate(returnInfo::ReReturnInfo) = first(returnInfo.state)

function returning(returnInfo::ReReturnInfo, restart::RestartState, verbose::Verbose)
  if atlevel(verbose, ReturnLevel)
  	println(returnInfo, restart)
  end
  reverse!(returnInfo)
  collect!(returnInfo)
  if monitorable(returnInfo.runInfo)
  	reverse!(returnInfo.runInfo)
  	collect!(returnInfo.runInfo)
  end
  returnInfo
end

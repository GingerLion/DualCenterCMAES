#---------------------------------------------
# abstract  System
#---------------------------------------------
# Interface fields
#  maxGen::Integer
#  maxEvals::Integer

#---------------------------------------------
# abstract  State
#---------------------------------------------
# Interface fields
#  gen::Integer
#  evalCount::Integer
#  popn::Population
#  status::Symbol
#  best::Tuple

status(state::State) = state.status
status_(state::State) = state.status_shadow
gen(state::State) = state.gen
gen_(state::State) = state.gen_shadow
currentgen(state::State) = state.gen
currentgen_(state::State) = state.gen_shadow
incgen!(state::State) = (state.gen += 1)
incgen!_(state::State) = (state.gen_shadow += 1)
decgen!(state::State) = (state.gen -= 1)
decgen!_(state::State) = (state.gen_shadow -= 1)
evals(state::State) = state.evalCount
evals_(state::State) = state.evalCount_shadow
localevals(state::State) = state.evalCount
localevals_(state::State) = state.evalCount_shadow
incevals!(state::State, evalsPerGen::Integer) = (state.evalCount = state.evalCount + evalsPerGen)
incevals!_(state::State, evalsPerGen::Integer) = (state.evalCount_shadow = state.evalCount_shadow + evalsPerGen)
evolvable(state::State) = (state.status == :evolve)
evolvable!(state::State) = (state.status = :evolve)
evolvable_(state::State) = (state.status_shadow == :evolve)
evolvable!_(state::State) = (state.status_shadow = :evolve)
completed(state::State, sys::System, f::Fitness) = found(state, f) || gtmaxgen(state, sys)
completed_(state::State, sys::System, f::Fitness) = found_(state, f) || gtmaxgen(state, sys)
gtmaxgen(state::State) = (state.gen > system(state).maxGen)
gtmaxgen!(state::State) = if gtmaxgen(state, system(state)) state.status = :max_gen end

gtmaxevals(evals::Int, sys::System) = (evals > sys.maxEvals)
gtmaxevals(state::State) = (evals(state) > system(state).maxEvals)
gtmaxevals_(state::State) = (evals_(state) > system(state).maxEvals)
gtmaxevals!(state::State) = if gtmaxevals(state) state.status = :max_evals; state.stopEvals = true end
gtmaxevals!_(state::State) = if gtmaxevals_(state) state.status_shadow = :max_evals; state.stopEvals_ = true end

function gtmaxevals!(state::State, restart::RestartState)
	if gtmaxevals(evals(restart), system(state))
		state.status = :max_evals
		state.stopEvals = true
	end
	if gtmaxevals(evals_(restart), system(state))
		state.status_shadow = :max_evals
		state.stopEvals_ = true
	end
end

system(state::State) = state.sys
population(state::State) = state.popn
population_(state::State) = state.sOffspring_shadow
popnsize(state::State) = popnsize(population(state))
popnsize_(state::State) = popnsize(population_(state))
chrlength(state::State) = chrlength(population(state))
chrlength_(state::State) = chrlength(population_(state))
dimensions(state::State) = chrlength(population(state))
dimensions_(state::State) = chrlength(population_(state))
maximizing(state::State) = maximizing(population(state))
maximizing_(state::State) = maximizing(population_(state))
direction(state::State) = direction(population(state))
direction_(state::State) = direction(population_(state))
best!(state::State) = (state.best = best(population(state)))
best(state::State) = state.best
best_(state::State) = state.best_shadow
best_overall(state::State) = state.best_overall
best_overall_(state::State) = state.best_overall_
bestchromosome(state::State) = bestchromosome(best(state))
bestchromosome_(state::State) = bestchromosome(best_(state))
bestchroverall(state::State) = bestchromosome(best_overall(state))
bestchroverall_(state::State) = bestchromosome(best_overall_(state))
bestfitness(state::State) = bestfitness(best(state))
bestfitness_(state::State) = bestfitness(best_(state))
bestfitoverall(state::State) = bestfitness(best_overall(state))
bestfitoverall_(state::State) = bestfitness(best_overall_(state))

better(state1::State, state2::Union{State,Tuple}) = (maximizing(state1) ? bestfitness(state1) > bestfitness(state2)
																 			  : bestfitness(state1) < bestfitness(state2))
better_(state1::State, state2::Union{State,Tuple}) = (maximizing_(state1) ? bestfitness_(state1) > bestfitness(state2)
																			  : bestfitness_(state1) < bestfitness(state2))
better_overall(state1::State, state2::Union{State, Tuple}) = (maximizing(state1) ? bestfitoverall(state1) > bestfitness(state2)
																				   : bestfitoverall(state1) < bestfitness(state2))
better_overall_(state1::State, state2::Union{State, Tuple}) = (maximizing_(state1) ? bestfitoverall_(state1) > bestfitness(state2)
																					: bestfitoverall_(state1) < bestfitness(state2))
better(member::Tuple, state::State) = (maximizing(state) ? bestfitness(member) > bestfitness(state)
															: bestfitness(member) < bestfitness(state))


worse(state1::State, state2::Union{State,Tuple}) = (maximizing(state1) ? bestfitness(state1) < bestfitness(state2)
																 			  : bestfitness(state1) > bestfitness(state2))
worse(member::Tuple, state::State) = (maximizing(state) ? bestfitness(member) < bestfitness(state)
															: bestfitness(member) > bestfitness(state))
same(state1::State, state2::Union{State,Tuple}) = bestfitness(state1) == bestfitness(state2)
same(member::Tuple, state::State) = bestfitness(member) == bestfitness(state)

function found(state::State, fit::Fitness)
	if found_(state)
	 	abs(bestfitoverall(state) - fit.optimalValue) <= fit.epsilon && better_overall(state, best_overall_(state))
	else
		abs(bestfitoverall(state) - fit.optimalValue) <= fit.epsilon
	end
end
function found_(state::State, fit::Fitness)
	if found(state)
		abs(bestfitoverall_(state) - fit.optimalValue) <= fit.epsilon && better_overall_(state, best_overall(state))
	else
		abs(bestfitoverall_(state) - fit.optimalValue) <= fit.epsilon
	end
end

function hastocatchup(state::State)
	if found(state) && (status_(state) == :stop || status_(state) == :evolve) && better_overall(state, best_overall_(state))
		return :dualcenter
	elseif found_(state) && (status(state) == :stop || status(state) == :evolve) && better_overall_(state, best_overall(state))
		return :normal
	else
		return :none
	end
end

#decides which system shouuld continue to evolve when there is no found and a max_evals
function maxContinueSystem(state::State)
	if (status(state) == :stop && status_(state) == :max_evals)
		return :normal
	elseif (status(state) == :max_evals && status_(state) == :stop)
		return :dualcenter
	else
		return :none
	end
end
#returns true if one system is max_evals and one is stop
isMaxAndStop(state::State) = ((status(state) == :stop && status_(state) == :max_evals) || (status(state) == :max_evals && status_(state) == :stop)) ? true : false

found(state::State) = (state.status == :found)
found_(state::State) = (state.status_shadow == :found)

function found!(state::State, fit::Fitness)
 	if found(state)
 		found(state)
 	elseif found(state, fit)
		println("Normal found solution at gen = $(currentgen(state))")
 		state.status = :found
		state.stopEvals = true
	end
end

function found!_(state::State, fit::Fitness)
 	if found_(state)
 		found_(state)
 	elseif found_(state, fit)
		println("Shadow found solution at gen = $(currentgen_(state))")
		state.status_shadow = :found
		state.stopEvals_ = true
	end
end

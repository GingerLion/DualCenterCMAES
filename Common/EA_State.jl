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
gtmaxevals!(state::State) = if gtmaxevals(state) state.status = :max_evals end

function gtmaxevals!(state::State, restart::RestartState)
	if gtmaxevals(evals(restart), system(state))
		state.status = :max_evals
	end
end

system(state::State) = state.sys
population(state::State) = state.popn
popnsize(state::State) = popnsize(population(state))
chrlength(state::State) = chrlength(population(state))
dimensions(state::State) = chrlength(population(state))
maximizing(state::State) = maximizing(population(state))
direction(state::State) = direction(population(state))
best!(state::State) = (state.best = best(population(state)))
best(state::State) = state.best
best_(state::State) = state.best_shadow
bestchromosome(state::State) = bestchromosome(best(state))
bestchromosome_(state::State) = bestchromosome(best_(state))
bestfitness(state::State) = bestfitness(best(state))
bestfitness_(state::State) = bestfitness(best_(state))
better(state1::State, state2::Union{State,Tuple}) = (maximizing(state1) ? bestfitness(state1) > bestfitness(state2)
																 			  : bestfitness(state1) < bestfitness(state2))
better(member::Tuple, state::State) = (maximizing(state) ? bestfitness(member) > bestfitness(state)
															: bestfitness(member) < bestfitness(state))


worse(state1::State, state2::Union{State,Tuple}) = (maximizing(state1) ? bestfitness(state1) < bestfitness(state2)
																 			  : bestfitness(state1) > bestfitness(state2))
worse(member::Tuple, state::State) = (maximizing(state) ? bestfitness(member) < bestfitness(state)
															: bestfitness(member) > bestfitness(state))
same(state1::State, state2::Union{State,Tuple}) = bestfitness(state1) == bestfitness(state2)
same(member::Tuple, state::State) = bestfitness(member) == bestfitness(state)

found(state::State, fit::Fitness) = abs(bestfitness(state) - fit.optimalValue) <= fit.epsilon
found_(state::State, fit::Fitness) = abs(bestfitness_(state) - fit.optimalValue) <= fit.epsilon
found(state::State) = (state.status == :found)
found_(state::State) = (state.status_shadow == :found)

function found!(state::State, fit::Fitness)
 	if found(state)
 		found(state)
 	elseif found(state, fit)
 		state.status = :found
	end
end

function found!_(state::State, fit::Fitness)
 	if found_(state)
		println("Shadow found solution first!")
 		found_(state)
 	elseif found_(state, fit)
		println("Shadow found solution first!")
		state.status_shadow = :found
		#state.status = :found
	end
end

#-----------------------------------------------------------------------------
# abstract type Population
#-----------------------------------------------------------------------------
# base instance variables - assumed to exist in all subclasses
#   members
#   fitness
#   direction::Symbol
#-----------------------------------------------------------------------------

# Interface:
#   chrlength, popnsize, 
#   members, members_, fitness, fitness_, direction,
#   minimizing, maximizing, minimize (depreciated), maximize (depreciated)
#   maximum, minimum, indbest, best, bestchromosome, bestfitness
#   getindex, setindex!, mcat, pcat, 
#   +, *,


chrlength(popn::Population) = size(popn.members,1)
popnsize(popn::Population) = size(popn.members,2)
direction(popn::Population) = popn.direction
members(nonPopn) = nonPopn
members(popn::Population) = copy(popn.members)
members_(popn::Population) = popn.members
fitness_(popn::Population) = popn.fitness
fitness(popn::Population) = copy(popn.fitness)
maximum(popn::Population) = popn[:both, argmax(popn)]
minimum(popn::Population) = popn[:both, argmin(popn)]
minimize(popn::Population) = (popn.direction == :min)
maximize(popn::Population) = (popn.direction == :max)
minimizing(popn::Population) = (popn.direction == :min)
maximizing(popn::Population) = (popn.direction == :max)
indbest(popn::Population) = minimizing(popn) ? argmin(popn) : argmax(popn)
best(popn::Population) = minimize(popn) ? minimum(popn) : maximum(popn)
bestchromosome(b::Tuple) = b[1]
bestfitness(b::Tuple) = b[2]
bestchromosome(popn::Population) = bestchromosome(best(popn))
bestfitness(popn::Population) = bestfitness(best(popn))
# setindex!(popn::Population, value, i) =  (popn.members[:,i] = value)      # should depreciate - shouldn't exist in general
+(values::Array, popn::Population) = values + members(popn)
+(popn::Population, values::Array) = members(popn) + values
*(values::Array, popn::Population) = values * members(popn)
*(popn::Population, values::Array) = members(popn) * values
*(values::Array, popn::Population) = values * members(popn)

function getindex(popn::Population, choice::Symbol, i)
  if (choice == :fit)
    popn.fitness[i]
  elseif (choice == :chr)
    popn.members[:,i]
  elseif (choice == :both)
    (popn[:chr,i], popn[:fit,i])
  else
    error("First index should be :chr_, :chr, :fit_, :fit, :both_ or :both. Instead it was $(choice)")
  end
end

# mcat is short for member concatenation

function mcat(popn::Population, rest...)
  allPopns = (popn, rest...)
  allMembers = map((x)->members(x), allPopns)
  hcat(allMembers...)
end


# pcat is short for population concatenation
# note: when using pcat, make sure that all populations have been evaluated
#       if not the case use RegularPopulation(mcat(popn1, popn2)) and then evaluate
# note: returns RegularPopuation as it would ruin the sort

function pcat(popn::Population, rest...)
  allPopns = (popn, rest...)
  allMembers = map((x)->members(x), allPopns)
  allFitness = map((x)->fitness(x), allPopns)
  combinedMembers = hcat(allMembers...)
  combinedFitness = vcat(allFitness...)
  RegularPopulation(popn, combinedMembers, combinedFitness)
end

#-----------------------------------------------------------------------------
# mutable struct RegularPopulation <: Population
#-----------------------------------------------------------------------------
#   members
#   fitness
#   direction::Symbol
#   evaluated::Bool
#-----------------------------------------------------------------------------

###  Constructors

RegularPopulation(members; direction = :min) = RegularPopulation(members, Array{Real}(undef, popnsize(members)), direction, false)
RegularPopulation(members, fitness; direction = :min) = RegularPopulation(members, fitness, direction, true)

RegularPopulation(m::Member; direction = :min) = RegularPopulation(chromosome(m), fitness(m), direction, true)

function RegularPopulation(similarPopn::Population, m::Member; direction = :min)
  RegularPopulation(chromosome(m), fitness(m), similarPopn.direction, true)
end

#called from Noise.jl      (n.values, members)
function RegularPopulation(similarPopn::Population, members)
  RegularPopulation(members, Array{Real}(undef, popnsize(members)), similarPopn.direction, false)
end
#called from Noise.jl       (n.values, members, fitness)
function RegularPopulation(similarPopn::Population, members, fitness)
  RegularPopulation(members, fitness, similarPopn.direction, true)
end

function RegularPopulation(chrLength::Integer, popnSize::Integer; direction = :min)
  RegularPopulation(Array{Real}(undef, chrLength, popnSize), Array{Real}(undef, popnSize), direction, false)
end

# function RegularPopulation(members, objfn::Function; direction = :min)
#   popn = RegularPopulation(members, direction = direction)
#   evaluate!(popn, objfn)
#   popn
# end

function RegularPopulation(members, fit::Fitness; direction = :min)
  popn = RegularPopulation(members, direction = direction)
  evaluate!(popn, fit)
  popn
end

function RegularPopulation(individual::Vector, popnSize::Integer; direction = :min)
  chrLength = length(individual)
  locus = individual[1]
  members = Array{typeof(locus), 2}(undef, chrLength, popnSize)
  for i = 1:popnSize
    members[:,i] = individual
  end
  RegularPopulation(members, direction = direction)
end

function RegularPopulation(individual::Vector, popnSize::Integer, objFn::Function; direction = :min)
  chrLength = length(individual)
  locus = individual[1]
  members = Array{typeof(locus), 2}(undef, chrLength, popnSize)
  fit = objFn(individual)
  for i = 1:popnSize
    members[:,i] = individual
  end
  RegularPopulation(members, fill(fit, popnSize); direction = direction)
end

# This is not the same as copy() or deepcopy()
# While kept in sorted order, the truncated population size becomes fixed
function RegularPopulation(popn::SortedPopulation)
  RegularPopulation(members(popn), fitness(popn), direction = popn.direction)
end

####  Other methods that extend or override Population methods

function getindex(popn::RegularPopulation, i)
  if (evaluated(popn))
    RegularPopulation(popn[:chr, i], popn[:fit, i], direction(popn), true)
  else
    RegularPopulation(popn[:chr, i]; direction = direction(popn))
  end
end

# mcat! is not defined for Population, because should not be done on SortedPopulation
#   - would ruin the sort

mcat!(popn::RegularPopulation, members, rest...) =  (popn.member = mcat(popn, members, rest))

# evaluate is only defined for RegularPopulation because sortedPopulation must already be evaluated

function evaluate(popn::RegularPopulation, f::Fitness)
  fitness = Array{returntype(f)}(undef, popnsize(popn))
  if any(j -> (j == -Inf), popn[:chr,:])
      println("WARNING: found -Infs and making them NaN!!!")
      for i = 1:popnsize(popn)
        fitness[i] = NaN
      end
  else
      for i = 1:popnsize(popn)
        fitness[i] = f.objFn(popn[:chr, i])
      end
  end
  fitness
end

function evaluate!(popn::RegularPopulation, f::Fitness)
  popn.evaluated = true
  popn.fitness = evaluate(popn, f)
end

argmax(popn::RegularPopulation) = argmax(popn.fitness)
argmin(popn::RegularPopulation) = argmin(popn.fitness)
optimize_as(popn::RegularPopulation, direction::Symbol) = (popn.direction == direction)
evaluated(popn::RegularPopulation) = popn.evaluated
evaluated!(popn::RegularPopulation, value = true) = (popn.evaluated = value)


function setindex!(popn::RegularPopulation, value, choice::Symbol, i::Integer)
  if (choice == :fit)
    popn.fitness[i] = value
  elseif (choice == :chr)
    popn.evaluated = false
    popn.members[:,i] = value
  elseif (choice == :both)
    (popn.members[:,i], popn.fitness[i]) = value
  else
    error("First index should be :chr, :fit, or :both. Instead it was $(choice)")
  end
end

setindex!(popn::RegularPopulation, value, i) =  (popn[:chr, i] = value)
setindex!(popn::RegularPopulation, value::Tuple, i) =  (popn[:both, i] = value)

# note: currently only applied to RegularPopulations as + or *  would change the fitness so ruin the sort
+(popn1::RegularPopulation, popn2::Population) = RegularPopulation(popn1, popn1.members + popn2.members)
*(popn1::RegularPopulation, popn2::Population) = RegularPopulation(popn1, popn1.members * popn2.members)

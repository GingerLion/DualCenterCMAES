#-----------------------------------------------------------------------------
# struct SortStructure      # used by other populations to match sort order of the SortedPopulation
#-----------------------------------------------------------------------------
#   fitness
#   sortOrder::Vector       # sort for members [1:λ], used for reseting μ in a truncation where μ_new > μ_old
#   index::Vector           # sort for members [1:μ], used for indexing into the members
#                           #   -  naturally causes truncated out of bound errors
#   lt::Function
#   direction::Symbol
#-----------------------------------------------------------------------------
#   used in conjunction with congruously structred RegularPopulations and the SortedPopulation constructors
#      to "sort" the congruous population in the same sort order with only a constant overhead
#-----------------------------------------------------------------------------

function SortStructure(popn::SortedPopulation)
  SortStructure(copy(popn.fitness),
                copy(popn.sortOrder),
                copy(popn.index),
                popn.lt,
                popn.direction)
end

#-----------------------------------------------------------------------------
# mutable struct SortedPopulation <: Population
#-----------------------------------------------------------------------------
#   members                 # holds the members from the original population - should not be modified
#   fitness
#   sortOrder::Vector       # sort for members [1:λ], used for reseting index in a truncation where μ_new > μ_old
#   index::Vector           # sort for members [1:μ], used for indexing into the members
#   lt::Function
#   direction::Symbol
#-----------------------------------------------------------------------------

function SortedPopulation(popn::RegularPopulation, lt::Function)
  index = sortperm(popn.fitness; lt = lt)
  SortedPopulation(popn.members, popn.fitness, copy(index), copy(index), lt, popn.direction)
end

SortedPopulation(popn::RegularPopulation) = SortedPopulation(popn, minimize(popn) ? (<) : (>))

function SortedPopulation(members, structure::SortStructure)
  SortedPopulation(members, structure.fitness, copy(structure.sortOrder),
                   copy(structure.index), structure.lt, structure.direction)
end

function SortedPopulation(popn::RegularPopulation, structure::SortStructure)
  SortedPopulation(popn.members, structure)
end

evaluated(popn::SortedPopulation) = true

function truncate!(popn::SortedPopulation, μ)
  @assert μ <= maxsize(popn) "truncated population size must be no larger then originating population"
  popn.index = popn.sortOrder[1:μ]
  #println("sort order = $(popn.index)\n")
  SortStructure(popn) #update SortStructure
end

popnsize(popn::SortedPopulation) = length(popn.index)
indmax(popn::SortedPopulation) = 1
indmin(popn::SortedPopulation) = 1
maxsize(popn::SortedPopulation) =  length(popn.sortOrder)
evaluated(popn::SortedPopulation) = true
evaluated!(popn::SortedPopulation, value = true) = (value == true) ? true : error("Sorted populations are already evaluated - they cannot be set to 'evaluation = false'")
fitness(popn::SortedPopulation) =  popn[:fit,:]
members(popn::SortedPopulation) =  popn[:chr,:]
#getindex(popn::SortedPopulation, i) = popn[:chr, i]

function getindex(popn::SortedPopulation, i)
  fit = popn[:fit, i]
  if !(typeof(fit)<:Vector) fit = [fit] end #fit = [fit] changes fit from Float64 to Array{Float64,1} with 1 elment
  SortedPopulation(popn[:chr, i], popn[:fit, i], collect(1:length(i)), collect(1:length(i)), popn.lt, direction(popn))
end

function getindex(popn::SortedPopulation, choice::Symbol, i)
  if (choice == :fit)
    popn.fitness[popn.index[i]]
  elseif (choice == :chr)
    popn.members[:,popn.index[i]]
  elseif (choice == :both)
    (popn[:chr,i], popn[:fit,i])
  else
    error("First index should be :chr, :fit, or :both. Instead it was $(choice)")
  end
end

best(popn::SortedPopulation) = popn[:both, 1]

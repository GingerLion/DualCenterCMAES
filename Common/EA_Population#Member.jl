#-----------------------------------------------------------------------------
# struct Member
#-----------------------------------------------------------------------------
### created from the population, not for the population,
###    which is based on a 2d matrix of real values
# struct Member
#   chromosome
#   fitness
# end

Member(pair::Tuple) =  Member(pair...)
Member(chr, objfn::Function) = Member(chr, objfn(chr))
chromosome(m::Member) = m.chromosome
fitness(m::Member) = m.fitness
print(m::Member) = print("{$(fitness(m))}\t$(chromosome(m))")

#-----------------------------------------------------------------------------
# helper functions for the members field of a Population
#-----------------------------------------------------------------------------

chrlength(members::Array) = size(members,1)
popnsize(members::Array) = size(members,2)

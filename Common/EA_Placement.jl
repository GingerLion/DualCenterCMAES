function place!(popn::RegularPopulation, member::Member; locate = indworst)
	loc = locate(popn, member)
	if loc > 0
		popn[loc] = member
	end
	loc
end

function worse(popn::Population, member::Member; loc = 1)
  (maximizing(popn) ? popn[:fit, loc] < fitness(member) 
            : popn[:fit, loc] > fitness(member))
end

indworst(popn::RegularPopulation) = maximizing(popn) ? indmin(popn) : indmax(popn)
indworst(popn::SortedPopulation) = popnsize(popn)

function indworst(popn::Population, member::Member)
	loc = indworst(popn)
	if worse(popn, member; loc = loc)
		loc
	else
		0
	end
end

replaceworst!(popn::RegularPopulation, member::Member) = place!(popn, member)

function replaceworst!(popn::SortedPopulation, member::Member)
	loc = indworst(popn)
	if worse(popn, member; loc = loc)
		internalLoc = popn.index[loc]
		popn.members[:, loc] = chromosome(member)
		popn.fitness[loc] = fitness(member)
		internalLoc
	else 
		0
	end
end

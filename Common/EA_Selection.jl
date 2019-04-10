randompermute(range) = sample(range, length(range); replace = false)

let extrasCount = 0, totalExtras = 0, fCurrentCount = 0, bCurrentCount = 0, fullCount = 0
  global function extra!(eCount, startCount)
    fCurrentCount = 1
    fullCount = bCurrentCount = startCount
    totalExtras = extrasCount = eCount
  end

  global function nextextra(loaded = :front)
    (loaded == :front) ? nextextra_front() : nextextra_back()
  end

  global function nextextra_front()
  	extra = (fCurrentCount > totalExtras) ? 0 : 1
    fCurrentCount += 1
    if fCurrentCount > fullCount
    	fCurrentCount = 1
    end
    extra
  end

  global function nextextra_back()
  	extra = (bCurrentCount > totalExtras) ? 0 : 1
    bCurrentCount -= 1
    if bCurrentCount == 0
    	bCurrentCount = fullCount
    end
    extra
  end
end

let segmentSize = 2
  global function segmentsize!(λ::Int, μ::Int)
    (segmentSize, extrasCount) = divrem(λ, μ)
    extra!(extrasCount, μ)
    segmentSize
  end

  global function segmentsize(loaded = :front)
    segmentSize + nextextra(loaded)
  end
end

function indexorder(μ::Int, beginBinding::Symbol)
  if (beginBinding == :leastFit)
  	indexOrder = μ:-1:1
  elseif (beginBinding == :mostFit)
  	indexOrder = 1:μ
  else
  	indexOrder = randompermute(1:μ)
  end
  indexOrder
end

function interleave(parents::Population, offspring::Population;
	                beginBinding = :mostFit, loaded = :mostFit)
  μ = popnsize(parents)
  λ = popnsize(offspring)
  @assert μ <= λ "the parents population must not have more members than the offspring population"
  segmentsize!(λ, μ)
  segmentSize = Array{Int}(μ)
  members = Array{Real}(chrlength(parents), λ + μ)
  fitness = Array{Real}(λ + μ)
  parentsOrder = Array{Int}(μ)
  offspringOrder = Array{Int}(λ)

  distances = pairwise(Euclidean(), members_(parents), members_(offspring))   # distances will be size(μ,λ)

  c = 1; p = 1; i = 1

  for parent in indexorder(μ, beginBinding)
    segSize = segmentSize[p] = segmentsize(loaded)
  	parentsOrder[p] = parent;   (p += 1)
  	(members[:,i], fitness[i]) = parents[:both, parent];   (i += 1)
  	for j = 1:segSize
    	child = indmin(distances[parent,:])
    	offspringOrder[c] = child;   (c += 1)
    	distances[:, child] = fill(Inf, μ)   # set dist to Inf so this child will not be chosen again
  		(members[:,i], fitness[i]) = offspring[:both, child];   (i += 1)
	  end
  end
  (RegularPopulation(parents, members, fitness), segmentSize, parentsOrder, offspringOrder)
end

function matchedinterleave(parents::Population, offspring::Population, segmentSize, pOrder, cOrder)
  if evaluated(parents) && evaluated(offspring)
    matchedinterleave_e(parents, offspring, segmentSize, pOrder, cOrder)
  else
    matchedinterleave_n(parents, offspring, segmentSize, pOrder, cOrder)
  end
end

function matchedinterleave_e(parents::Population, offspring::Population, segmentSize, pOrder, cOrder)
  members = Array{Real}(chrlength(parents), popnsize(parents) + popnsize(offspring))
  fitness = Array{Real}(chrlength(parents), popnsize(parents) + popnsize(offspring))
  i = 1; c = 1

  for p = 1:popnsize(parents)
    (members[:,i], fitness[i]) = parents[:both, pOrder[p]]
    i += 1
    for j = 1:segmentSize[p]
      (members[:,i], fitness[i]) = offspring[:chr, cOrder[c]]#bug?
      i += 1; c += 1
    end
  end

  RegularPopulation(parents, members, fitness)
end

function matchedinterleave_n(parents::Population, offspring::Population, offspringPerSeg, pOrder, cOrder)
  members = Array{Real}(chrlength(parents), popnsize(parents) + popnsize(offspring))
  i = 1; c = 1

  #println("sizes: members = $(size(members)), parents = $(popnsize(parents)), offspring = $(popnsize(offspring))")
  #println("cOrder = $cOrder, pOrder = $pOrder, segSize = $offspringPerSeg")
  for p = 1:popnsize(parents)
    members[:,i] = parents[:chr, pOrder[p]]
    i += 1
    for j = 1:offspringPerSeg[p]
      members[:,i] = offspring[:chr, cOrder[c]]
      i += 1; c += 1
    end
  end
  RegularPopulation(parents, members)
end

function truncation(popn::RegularPopulation, slotSize::Array, truncationSize::Int, segCount::Integer)
  members = Array{Real}(chrlength(popn), truncationSize)
  fitness = Array{Real}(truncationSize)
  sortOrder = Array{Int}(truncationSize)

  segmentsize!(truncationSize, segCount)
  segEnd = mEnd = 0

  for i = 1:segCount
    segSize = segmentsize()
    segStart = segEnd + 1
    segEnd += segSize

    membersCount = sum(map((x)->slotSize[x], segStart:segEnd))
    mStart = mEnd + 1
    mEnd += membersCount
	  sPopn = SortedPopulation(popn[mStart:mEnd])
	  sortOrder[segStart:segEnd] = (sPopn.sortOrder .+ mStart)[1:segSize]
    (members[:,segStart:segEnd], fitness[segStart:segEnd]) = sPopn[:both, 1:segSize]
  end

  (RegularPopulation(popn, members, fitness), sortOrder)
end

function truncation(popn::RegularPopulation, truncationSize::Int, segCount::Int)
  members = Array{Real}(chrlength(popn), truncationSize)
  fitness = Array{Real}(truncationSize)
  sortOrder = Array{Int}(truncationSize)

  segmentsize!(truncationSize, segCount)
  segEnd = 0

  for i = 1:segCount
    segSize = segmentsize()
    segStart = segEnd + 1
    segEnd += segSize

	sPopn = SortedPopulation(popn[segStart:segEnd])
	sortOrder[segStart:segEnd] = (sPopn.sortOrder .+ segStart)[1:segSize]
    (members[:,segStart:segEnd], fitness[segStart:segEnd]) = sPopn[:both, 1:segSize]
  end

  (RegularPopulation(popn, members, fitness), sortOrder)
end

function truncation(popn::RegularPopulation, truncationSize::Int)
  sPopn = SortedPopulation(popn)
  truncate!(sPopn, truncationSize)

  (sPopn, SortStructure(sPopn))
end

function truncation(sPopn::SortedPopulation, truncationSize::Int)
  sPopn = deepcopy(sPopn)
  truncate!(sPopn, truncationSize)

  (sPopn, SortStructure(sPopn))
end

function matchedtruncation(popn::RegularPopulation, order::Array)
  if evaluated(popn)
    matchedtruncation_e(popn, order)
  else
    matchedtruncation_n(popn, order)
  end
end

function matchedtruncation_e(popn::RegularPopulation, order::Array)
  members = Array{Real}(chrlength(popn), length(order))
  fitness = Array{Real}(length(order))

  for i = 1:length(order)
  	members[:, i] = popn[:chr, order[i]]
    fitness[i] = popn[:fit, order[i]]
  end

  RegularPopulation(popn, members, fitness)
end

function matchedtruncation_n(popn::RegularPopulation, order::Array)
  members = Array{Real}(chrlength(popn), length(order))

  for i = 1:length(order)
    members[:, i] = popn[:chr, order[i]]
  end

  RegularPopulation(popn, members, fitness)
end

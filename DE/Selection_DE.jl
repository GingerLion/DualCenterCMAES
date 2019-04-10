##########################################################################################################
###
###   Selection
###
##########################################################################################################

function slotselection(target::Population, samples::Population)
  allSamples = pcat(target, samples)
  best(allSamples)
end

slotselection(popn::Population, samples::Array, de::DE_System) = slotselection(popn, samples)


function slotselection(popn::Population, samples::Array)
  members = Array{Real, 2}(undef, chrlength(popn), popnsize(popn))
  fitness = Array{Real, 1}(undef, popnsize(popn))
  for i = 1:popnsize(popn)
    (members[:, i], fitness[i]) = slotselection(popn[i], samples[i])
  end
  RegularPopulation(members, fitness; direction = popn.direction)
end

function combine(samples::Array)
  chrlength = chrlength(samples[1])
  sampleCount = sum(map((s) -> popnsize(s), samples))
  allSamples = Array{Real, 2}(undef, chrlength, sampleCount)
  allFitness = Array{Real, 1}(undef, sampleCount)
  j = 1
  for s in samples
    for i = 1:popnsize(samples[i]) 
      (allSamples[:,j], allFitness[j]) = s[:both, i]
      j += 1
    end
  end
  RegularPopulation(allSamples, allFitness)
end

function interleave(popn::Population, samples::Array; indices = 1:popnsize(popn))
  chrLength = chrlength(popn)
  popnSize = popnsize(popn)
  sampleCount = sum(map((s) -> popnsize(s), samples))
  allMembers = Array{Real, 2}(undef, chrLength, sampleCount + popnSize)
  allFitness = Array{Real, 1}(undef, sampleCount + popnSize)
  slotSize = Array{Integer, 1}(undef, popnSize)
  k = 1
  for i = 1:popnSize
    slotSize[indices[i]] = popnsize(samples[indices[i]])
    (allMembers[:,k], allFitness[k]) = popn[:both, indices[i]]
    k += 1
    for j = 1:slotSize[indices[i]] 
      (allMembers[:,k], allFitness[k]) = samples[indices[i]][:both, j]
      k += 1
    end
  end
  interleavedPopn = RegularPopulation(allMembers, allFitness; direction = popn.direction)
  (interleavedPopn, slotSize .+ 1)
end

# partial assertion:  groupCount | popnsize(popn)
#     if popnsize(popn) mod groupCount > 0, add the extra "slots" to other groups
#     so groupSize = floor(popnsize(popn) / groupCount) or ceil(popnsize(popn) / groupCount)
#     and the number of groups equals groupCount exactly
# assertion: sum(slotsize) == length(newPopn)

randompermute(range) = sample(range, length(range); replace = false)

function truncationselection(popn::Population, samples::Array, de::DE_System)
  if de.segCount == 1
    truncationselection(popn, samples)
  elseif de.segCount == popnsize(popn)
    slotselection(popn, samples)
  elseif 1 < de.segCount < popnsize(popn)
    truncationselection(popn, samples, de.segCount; randomizeSlots = de.randomizeSlots)
  else
    error("For truncationselection(), 0 < groupCount <= $(popnSize) (popnSize). Instead groupCount = $(groupCount)")
  end
end

function truncationselection(popn::Population, samples::Array)
  (newPopn, slotSize) = interleave(popn, samples)
  sortedPopn = SortedPopulation(newPopn)
  truncate!(sortedPopn, popnsize(popn))
  sortedPopn
end

function truncationselection(popn::Population, samples::Array, segCount::Integer; randomizeSlots = false)
  chrLength = chrlength(popn)
  popnSize = popnsize(popn)

  members = Array{Real}(chrLength, popnSize)
  fitness = Array{Real}(popnSize)

  segSizeBase = floor(Int, popnSize/segCount)
  extrasCount = popnSize % segCount

  indices = randomizeSlots ? randompermute(1:popnSize) : 1:popnSize
  (newPopn, slotSize) = interleave(popn, samples; indices = indices)

  segEnd = mEnd = 0
  
  for i = 1:segCount
    extra = (extrasCount > 0) ? 1 : 0
    extrasCount -= extra

    segSize = segSizeBase + extra
    segStart = segEnd + 1
    segEnd += segSize

    membersCount = sum(map((x)->slotSize[x], segStart:segEnd))
    mStart = mEnd + 1
    mEnd += membersCount

    (members[:,segStart:segEnd], fitness[segStart:segEnd]) = SortedPopulation(newPopn[mStart:mEnd])[:both, 1:segSize]
  end

  RegularPopulation(members, fitness; direction = popn.direction)
end


#-------------------------------------
# Choosing Sampling-Selection function
#-------------------------------------
function select_fn(popnSize, sel, slope)
  if sel == :uniform
    sel_fn = runiform_fn(popnSize)
  elseif sel == :linear
    sel_fn = ((slope == 2) ? rlinear_fn(popnSize)
                           : rlinear_fn(popnSize; slope = slope))
  else
    error("sampling selection must be :uniform or :linear. Instead $sel used.")
  end
  sel_fn
end


#-------------------------------------
# Random Sampling-Selection
#-------------------------------------
runiform_fn(popnSize::Int) = () -> rand(1:popnSize)


#-------------------------------------
# Rank Sampling-Selection
#-------------------------------------

# produces linear rank selection positions, but is not linked into selection mechanism
# 0 <= slope <= 1
function rlinear_fn(popnSize::Int; slope = (popnSize - 1.0)/(popnSize + 1.0))
  n = popnSize
  m = -2.0/n/(n-1.0)*slope
  a = m * n
  b = (2 - m*n^2)

  function ()
    c = -2 * popnSize * rand()
    ceil(Int, minimum(roots(Poly([c,b,a]))))
  end
end

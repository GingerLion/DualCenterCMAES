
###########################################################
# Sampling Methods

# struct MethodShare
#   share::Vector{Numbers}
#   procedure::Symbol
# end

# struct SamplingMethod
#   methods::Vector{Function}
#   share::Vector{MethodShare}
#   placement::Vector{Union{Integer, Symbols}}
# end

# smlist = samplingmethods(sm1, sm2, sm3, ...)   is the same as
# smlist = samplingmethods(SamplingMethod(sm1), SamplingMethod(sm2), SamplingMethod(sm3), ...)

# smlist = samplingmethods(SamplingMethod(methods = [sm1, sm2, sm3], MethodShare([popnSize - 1, 1]), :exact)),
#                          SamplingMethod(sm4), SamplingMethod(sm5), ...)
# smlist = samplingmethods(SamplingMethod([sm1, sm2, sm3], MethodShare([1/3, 1/2, 1/6]),
#                          SamplingMethod(sm4), SamplingMethod(sm5), ...)
# smlist = samplingmethods(SamplingMethod([sm1, sm2, sm3], [1, 1, 1, 2, 2, 3, 3, 3, 3, 2, 2, 1]), 
#                          SamplingMethod(sm4), SamplingMethod(sm5), ...)
# smlist = samplingmethods(SamplingMethod([sm1, sm2, sm3], [1, x, x, 2, x, x, 3, x, x, 1, 2, 3]), 
#                          SamplingMethod(sm4), SamplingMethod(sm5), ...)
# smlist = samplingmethods(SamplingMethod([sm1, sm2, sm3], 
#                                         MethodShare([1/3, 1/2, 1/6], :onAverage),
#                                         [1, x, x, 2, x, x, 3, x, x, 1, 2, 3]), 
#                          SamplingMethod(sm4), SamplingMethod(sm5), ...)
# smlist = samplingmethods(SamplingMethod(([sm1, sm2, sm3], [sm4, sm5]), 
#                                          (MethodShare([1/3, 1/2, 1/6], :onAverage), MethodShare([n - 1, 1], :exact)),
#                                          [1, :x1, :x1, 2, :x2, :x2, 3, :x1, :x2, 1, 2, 3]),
#                          SamplingMethod(sm4), SamplingMethod(sm5), ...)
#
# implement above at a later date

# Simpler version where the sampling methods are the same for each target,
# but begins to move to code where you could change the methods for different targets
# Creates an sm of type Array{Function, 2} instead of Array{Function} that is currently being used
#
# original use
# sm = samplingmethods(npop, classicsampling_uxover, classicsampling_uxover, directsampling)
#
# function samplingmethods(popnSize, sm1::Function, smx::Vararg{Function})
#   altSize = length(smx) + 1
#   sm_template = Array{Function}(altSize)
#   sm_template[1] = sm1
#   sm_template[2:altSize] = smx
#   sm = Array{Function}(altSize, popnSize)
#   for i = 1:popnSize
#     sm[:, i] = sm_template
#   end
# end

# struct SlotSamplers
#   samplers::Dict
#   function SlotSamplers()
# end

function SlotSamplers(popnSize::Int, sb::DE_SamplingParms)
  localSize = ceil(Int, popnSize * sb.localDonor)      
  donorselection = select_fn(localSize, sb.donorSel, sb.donorSlope)
  parentselection = select_fn(popnSize, sb.parentSel, sb.parentSlope)

  s = SlotSamplers()
  s[:classic] = classicsamplinguxover_fn(sb, donorselection, parentselection)
  s[:classic_direct] = classicsamplingdirect_fn(sb, donorselection, parentselection)
  s[:best] = bestsamplinguxover_fn(sb, parentselection)
  s[:best_direct] = bestsamplingdirect_fn(sb, parentselection)
  s[:direct] = directsampling_fn(sb, parentselection)
  s
end

samplertemplate(s::SlotSamplers, key::Symbol) = s[key]

function samplertemplate(s::SlotSamplers, keys::Array{Symbol})
  map(keys) do x
    s[x]
  end
end

getindex(s::SlotSamplers, key::Symbol) = s.samplers[key]

setindex!(s::SlotSamplers, value::Function, key::Symbol) = (s.samplers[key] = value)


# Classic mode

function ClassicParents(popnSize::Int, target_index::Int, donorselect::Function, parentselect::Function)
  parent1_index = parentselect()
  parent2_index = (parent1_index + parentselect() - 1) % popnSize + 1
  donor_index = (target_index + donorselect() - 1) % popnSize + 1
  ClassicParents(parent1_index, parent2_index, donor_index, target_index)
end

function createalternate(popn::Population, parents::ClassicParents, diffWeight::Float64)
    popn[:chr, parents.base] + diffWeight * (popn[:chr, parents.diffIndv1] - popn[:chr, parents.diffIndv2])
end

function createalternate(popn::Population, parents::BestParents, diffWeight::Float64)
    popn[:chr, parents.best] + diffWeight * (popn[:chr, parents.diffIndv1] - popn[:chr, parents.diffIndv2])
end

function createcombined(target::Vector, alternate::Vector, CR::Real)
  offspring = copy(alternate)
  chrLength = length(target)
  force = rand(1:chrLength)
  for m = 1:chrLength
    if rand() > CR && m != force
        offspring[m] = target[m]    # update matrix         (L4)
    end
  end
  offspring
end

function classicsamplinguxover_fn(de::DE_SamplingParms, donorsel::Function, parentsel::Function)
  function (popn::Population, target_index::Integer)
      parents = ClassicParents(popnsize(popn), target_index, donorsel, parentsel)
      alternate = createalternate(popn, parents, de.diffWeight)
      target = popn[:chr, target_index]
      combined = createcombined(target, alternate, de.CR)
      (DE_Individual(combined), DE_Individual(alternate))
    end
end

function classicsamplingdirect_fn(de::DE_SamplingParms, donorsel::Function, parentsel::Function)
  function (popn::Population, target_index::Integer)  
    parents = ClassicParents(popnsize(popn), target_index, donorsel, parentsel)
    alternate = createalternate(popn, parents, de.diffWeight)
    (DE_Individual(alternate), DE_Individual(alternate))
  end
end


# Best mode

function BestParents(popnSize::Int, target_index::Int, best_index::Int, parentselect::Function)
  parent1_index = parentselect()
  parent2_index = (parent1_index + parentselect() - 1) % popnSize + 1
  BestParents(parent1_index, parent2_index, best_index, target_index)
end

function bestsamplinguxover_fn(de::DE_SamplingParms, parentselect::Function)
  function (popn::Population, target_index::Integer)
      parents = BestParents(popnsize(popn), target_index, indbest(popn), parentselect)
      alternate = createalternate(popn, parents, de.diffWeight)
      target = popn[:chr, target_index]
      combined = createcombined(target, alternate, de.CR)
      (DE_Individual(combined), DE_Individual(alternate))
    end
end

function bestsamplingdirect_fn(de::DE_SamplingParms, parentselect::Function)
  function (popn::Population, target_index::Integer)  
    parents = BestParents(popnsize(popn), target_index, indbest(popn), parentselect)
    alternate = createalternate(popn, parents, de.diffWeight)
    (DE_Individual(alternate), DE_Individual(alternate))
  end
end

# Direct mode

function DirectParents(popnSize::Int, target_index::Int, parentselect::Function)
  parent1_index = parentselect()
  parent2_index = (parent1_index + parentselect() - 1) % popnSize + 1
  DirectParents(parent1_index, parent2_index, target_index)
end

function createalternate(popn::Population, parents::DirectParents, diffWeight::Real)
  popn[:chr, parents.target] + diffWeight * (popn[:chr, parents.diffIndv1] - popn[:chr, parents.diffIndv2])
end

function directsampling_fn(de::DE_SamplingParms, parentselect::Function)
  function (popn::Population, target_index::Integer)  
    parents = DirectParents(popnsize(popn), target_index, parentselect)
    alternate = createalternate(popn, parents, de.diffWeight)
    (DE_Individual(alternate), DE_Individual(alternate))
  end
end

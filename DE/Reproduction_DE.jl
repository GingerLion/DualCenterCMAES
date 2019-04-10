##########################################################################################################
###
###   Reproduction
###
##########################################################################################################

## general DE evolution
genes(chr::DE_Individual) = chr.genes
length(chr::DE_Individual) = length(genes(chr))
getindex(chr::DE_Individual, i) = genes(chr)[i]
chrlength(samples::Array{DE_Individual}) = length(samples[1])
size(samples::Array{DE_Individual}) = (chrlength(samples), length(samples))

function RegularPopulation(samples::Array{DE_Individual})
  members = Array{Real, 2}(undef, chrlength(samples), length(samples))
  for i = 1:length(samples)
    members[:,i] = genes(samples[i])
  end
  RegularPopulation(members)
end

function reproduce(popn::Population, target_index::Integer, samplers::Array{Function})
  sampleCount = length(samplers)
  samples = Array{DE_Individual}(undef, sampleCount)
  alternates = Array{DE_Individual}(undef, sampleCount)
  for i = 1:sampleCount
    samples[i], alternates[i] = samplers[i](popn, target_index)
  end
  samplePop = RegularPopulation(samples)
  altPop = RegularPopulation(alternates)
  (samplePop, altPop)
end

# this version of producesamples to be used when samplingmethods are expanded 
#  so different slots can have different methods
# function producesamples(popn::Population, target_index::Integer, de::DE_Parms)
#   samples = Array{DE_Individual}(length(sm))
#   alternates = Array{DE_Individual}(length(sm))
#   for i = 1:length(sm)
#     samples[i], alternates[i] = de.sm[i, target_index](popn, target_index, de)
#   end
#   (RegularPopulation(samples), RegularPopulation(alternates))
# end

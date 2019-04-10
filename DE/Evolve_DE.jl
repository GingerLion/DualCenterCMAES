
function evolvepopn!(state::DE_State, fit::DE_Fitness)
  de = system(state)
  (state.popn, state.contenders) = evolve(population(state), fit, de)
  best!(state)
  found!(state, fit)              # must be run after best!(state)
  incgen!(state)
  incevals!(state, de.evalsPerGen)
end

function Contenders(popn::Population)
  targets = deepcopy(popn)
  samples = Array{RegularPopulation}(undef, popnsize(popn))
  alternates = Array{RegularPopulation}(undef, popnsize(popn))
  Contenders(targets, samples, alternates)
end

function evolve(popn::Population, fit::DE_Fitness, de::DE_System)
  c = Contenders(popn)

  for i = 1:popnsize(popn)
    c.samples[i], c.alternates[i] = reproduce(popn, i, de.samplers)
    evaluate!(c.samples[i], fit)
  end

  newPopn = de.selection(c.targets, c.samples, de)
  evaluated!(newPopn)
  (newPopn, c)
end

# this version of evolve to be used when samplingmethods are expanded 
#  so different slots can have different methods - will need to be refactored
# function evolve(popn::Population, sm::Array{Function, 2}, de::DE_System)
#   newPopn = RegularPopulation(chrlength(popn), popnsize(popn))
#   c = Contenders(popn)

#   for i = 1:popnsize(popn)
#     c.samples[i], c.alternates[i] = reproduce(popn, de.sm, de)
#     evaluate!(samples[i], fit.objFn, fit.direction)   
#     newPopn[:chr, i], newPopn[:fit, i] = slotselection(c.targets[:popn, i], samples[i])
#   end
#   (newPopn, c)
# end



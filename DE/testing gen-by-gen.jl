function consolidate(popnArray::Array{RegularPopulation})
	popnSize = length(popnArray)
	chrLength = chrlength(popnArray[1])
	members = Array{Float64,2}(chrLength, popnSize)
	if evaluated(popnArray[1])
		fitness = Array{Float64,1}( popnSize)
		for i = 1:popnSize
			members[:,i] = popnArray[i][:chr, 1]
			fitness[i] = popnArray[i][:fit, 1]
		end
		RegularPopulation(popnArray[1], members, fitness)
	else
		for i = 1:popnSize
			members[:,i] = popnArray[i][:chr, 1]
		end
		RegularPopulation(popnArray[1], members)
	end
end


sphere_fn(n; ε = 1.0e-10) = DE_Fitness(sphere, n, :min, 0.0, ε, 5.12, -5.12)
f = sphere_fn(20)

function de_setup(fit; popn_size = 100, max_gen = 100_000, selection = slotselection, segCount = 1,
                   returnValue = :summary, template = [:classic, :classic, :classic_direct, :direct], 
                   donorSel = :uniform, parentSel = :uniform, CR = 0.5)
  deBase = DE_Base(fit, max_gen, popn_size; selection = selection, segCount = segCount, donarSel = donorSel, parentSel = parentSel, CR = CR)
  samplers = samplertemplate(SlotSamplers(deBase),template)
  de = DE_System(deBase, samplers)
  (de, deBase, samplers)
end

de, deBase, samplers = de_setup(f, popn_size = 5, max_gen = 200, template = [:classic],
								selection = slotselection, segCount = 1,
                    			donorSel = :uniform, parentSel = :uniform, CR = 0.5)
state = DE_State(de, f)

p0 = population(state)
p0[:chr, :]
p0[:fit, :]

evolve!(state, f, de)

p1 = population(state)
p1.members
p1.fitness

c1 = state.contenders
s1 = consolidate(c1.samples)
a1 = consolidate(c1.alternates)
t1 = c1.targets
s1.members
a1.members
t1.members
s1.fitness
t1.fitness

evolve!(state, f, de)

p2 = population(state)
p2.members
p2.fitness

c2 = state.contenders
s2 = consolidate(c2.samples)
a2 = consolidate(c2.alternates)
t2 = c2.targets
s2.members
a2.members
t2.members
s2.fitness
t2.fitness

evolve!(state, f, de)

p3 = population(state)
p3.members
p3.fitness

c3 = state.contenders
s3 = consolidate(c3.samples)
a3 = consolidate(c3.alternates)
t3 = c3.targets
s3.members
a3.members
t3.members
s3.fitness
t3.fitness

de, deBase, samplers = de_setup(f, popn_size = 100, max_gen = 1000, template = [:classic],
								selection = slotselection, segCount = 1,
                    			donorSel = :uniform, parentSel = :uniform, CR = 0.5)
state = DE_State(de, f)

@time runDE!(state, de, f; verbose = false, returnValue = :summary)

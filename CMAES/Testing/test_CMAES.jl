#test population
p1 = RegularPopulation([1.0,2.0,3.0,4.0,5.0,6.0,7.0],13)
p2 = RegularPopulation(mcat(p1,p1))
p3 = p1 + p1
p1[1]
p2[:chr,1]
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
rosenbrock() = typeof(rosenbrock([0,0]))
evaluate(p1, rosenbrock)
evaluate!(p1, rosenbrock)
p4 = pcat(p1,p1)
p4[:chr,18]
p4[:fit,18]
mbr = [1 2 3; 5 6 7; 7 8 9; 10 11 12; 14 15 16]'
simplefitfn(x) = abs(mean(x) - 8)
simplefitfn() = typeof(mean([1,2]))
p5 = RegularPopulation(mbr)
evaluate!(p5, simplefitfn)
s5 = SortedPopulation(p5)
s5[1]
s5[2]
s5[:fit,1]
s5[:fit,2]
s5[:fit,3]
s5[:chr,3]
srtOdr = truncate!(s5,4)
p6 = RegularPopulation([1 2 3 4 5])
s6 = SortedPopulation(p6, srtOdr)
s6[3]

# Begin to go through evolve! one line at a time
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
rosenbrock() = typeof(rosenbrock([0,0]))

#test initial CMAES_Model
N = 2
startCenter = zeros(N)
mparms =  Model_Parms(N)
gen = 0
model = CMAES_Model(mparms, startCenter; gen = gen)

gen = 1
#test Noise
nE = SphericalNoise(mparms.N, mparms.λ)
members(nE)
nW = ShapedNoise(nE, model)

#test CMAES_Model and Noise
nOffspring = model + nW

# continue with evolve!
evaluate!(nOffspring, rosenbrock)

(offspring, E, W) = (nOffspring, nE, nW)
sOffspring = SortedPopulation(offspring)
sortStructure = truncate!(sOffspring, mparms.μ)

sW = sort(W, sortStructure)
flatten!(sW, weights(model))

# testing model update in slow motion
m1 = deepcopy(model)
center!(m1, weightedavg(sW))
stepsize!(m1, weightedavg(sW))
h_σ!(m1, gen)
covarmatrix!(m1, sW)
eigendecomp!(m1)

gen = 2
nE2 = SphericalNoise(mparms.N, mparms.λ)
members(nE2)
nW2 = ShapedNoise(nE2, m1)
nOffspring2 = m1 + nW2
evaluate!(nOffspring2, rosenbrock)
(offspring, W) = (nOffspring2, nW2)
sOffspring = SortedPopulation(offspring)
sortStructure = truncate!(sOffspring, mparms.μ)
sW = sort(W, sortStructure)
flatten!(sW, weights(model))
m2 = deepcopy(m1)
center!(m2, weightedavg(sW))
stepsize!(m2, weightedavg(sW))
h_σ!(m2, gen)
covarmatrix!(m2, sW)
eigendecomp!(m2)

# continue test of CMAES_Model
nModel = deepcopy(model)
update!(nModel, sW)

#test CMAES
N = 2
c = CMAES(rosenbrock, N, startingCenter = zeros(N), maxGen = 1000, optimalValue = 0)
evolve!(c)
c = CMAES(rosenbrock, N, startingCenter = zeros(N), maxGen = 1000, optimalValue = 0)
runcmaes!(c, sys, f)
foundsolution(c)
maxgenreached(c)
completed(c)


n = 2; optimalValue = 0.0
rosenbrock(x::Array{Float64, 1}) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
#rosenbrock() = typeof(rosenbrock([0,0]))
f = RealFitness(rosenbrock, n, :min, optimalValue, 1.0e-10)

sys = CMAES_System(n, f; maxGen = 1000, includeCenter = false, elitism = false)
c = CMAES_State(sys, f; startingCenter = zeros(n))
s = Storage()
runcmaes!(c, sys, f; storage = s, returnValue = :summary)
println(evals(c))
print(s, :best_member)
print(s, :new_center)
print(s, :new_σ)
println(evals(c))


sys = CMAES_System(n, f; maxGen = 1000, includeCenter = false, elitism = true)
c = CMAES_State(sys, f; startingCenter = zeros(n))
s = Storage()
runcmaes!(c, sys, f; storage = s)
println(evals(c))
print(s, :best_member)
print(s, :new_center)
print(s, :new_σ)
println(evals(c))

sys = CMAES_System(n, f; maxGen = 1000, includeCenter = true, elitism = false)
c = CMAES_State(sys, f; startingCenter = zeros(n))
s = Storage()
runcmaes!(c, sys, f; storage = s)
println(evals(c))
print(s, :best_member)
print(s, :new_center)
print(s, :new_σ)
println(evals(c))

sys = CMAES_System(n, f; maxGen = 1000, includeCenter = true, elitism = true)
c = CMAES_State(sys, f; startingCenter = zeros(n))
s = Storage()
runcmaes!(c, sys, f; storage = s)
println(evals(c))
print(s, :best_member)
print(s, :new_center)
print(s, :new_σ)
println(evals(c))

sys = CMAES_System(n, f; maxGen = 1000, includeCenter = false, elitism = true, segmentCount = 2)
c = CMAES_State(sys, f; startingCenter = zeros(n))
s = Storage()
runcmaes!(c, sys, f; storage = s)
println(evals(c))
print(s, :best_member)
print(s, :new_center)
print(s, :new_σ)
println(evals(c))

n = 3
testFn = generatetests(n)
sys = CMAES_System(n, testFn[:schwefel]; maxGen = 1000, includeCenter = false, elitism = false)
c = CMAES_State(sys, testFn[:schwefel]; startingCenter = zeros(n))
s = Storage()
runcmaes!(c, sys, testFn[:schwefel]; storage = s)

let n = 3
	testFn = generatetests(n)
	f = testFn[:rastrign]
	startPoint = ones(n)

	sys = CMAES_System(n, f; maxGen = 2000, includeCenter = false, elitism = false)
	c = CMAES_State(sys, f; startingCenter = startPoint)
	s = Storage()
	runcmaes!(c, sys, f; storage = s)
end

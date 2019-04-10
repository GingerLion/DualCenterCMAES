n = 30
deg = 30.0
btestFn = generatetests(n, deg; ε = 1.0e-5)
btestFn = generatetests(n; ε = 1.0e-5)
f = btestFn[:rastrigin]
f = btestFn[:elliptical]
f = btestFn[:ackley]
f = btestFn[:levy]

# sys = CMAES_System(n, f; maxGen = 10000, includeCenter = false, elitism = true, segmentCount = 5)
# sys = CMAES_System(n, f; maxEvals = 1_000_000, includeCenter = false, elitism = false, λ = 320)

sys = CMAES_System(n, f; maxEvals = 1_000_000, includeCenter = true, elitism = true)
rsys = CMAES_Restart(; η = 2)

# for most
cmaes = runEA(sys, f; center_init = ones(n), σ_init = 1.0, verbose = ReturnLevel());
ipop = runEA(sys, rsys, f; center_init = ones(n), σ_init = 1.0, verbose = RestartLevel(), monitoring = true);

# for levy and rosenbrock
cmaes = runEA(sys, f; center_init = zeros(n), σ_init = 1.0, verbose = false);
ipop = runEA(sys, rsys, f; center_init = zeros(n), σ_init = 1.0, verbose = false);

# Write output to file
outputPath = "/Users/mwineberg/Documents/2 Work/Projects/CMAES_DE"
rInfo = ipop.runInfo;
write_selectionsource(sys, rInfo; path = outputPath)

# run one generation at a time
d = Array{CMAES_State}(5)
d[1] = deepcopy(c)
for i = 2:5
	evolve!(c, sys, f)
	d[i] = deepcopy(c)
end

for i=1:5 
	println("d[$i].nW = $(d[i].nW)") 
end

for i=1:5 
	println("d[$i].sW = $(d[i].sW)") 
end

for i=1:5 
	println("d[$i].nOffspring = $(d[i].nOffspring)") 
end

for i=1:5 
	println("d[$i].sOffspring = $(d[i].sOffspring)") 
end

for i=1:5 
	println("d[$i].nModel.σ = $(d[i].nModel.σ)") 
end

for i=1:5 
	println("d[$i].nM.center = $(d[i].nModel.center)") 
end

for i=1:5 
	println("d[$i].pM.center = $(d[i].pModel.center)") 
end

for i=1:5 
	println("d[$i].pM.C = $(d[i].pModel.C)") 
end

for i=1:5 
	println("d[$i].pM.B = $(d[i].pModel.B)") 
end

for i=1:5 
	println("d[$i].pM.D = $(d[i].pModel.D)") 
end

for i=1:5 
	println("d[$i].pM.p_σ = $(d[i].pModel.p_σ)") 
end

for i=1:5 
	println("d[$i].pM.σ = $(d[i].pModel.σ)") 
end

for i=1:5 
	println("d[$i].nW = $(d[i].nW)") 
end

for i=1:5 
	println("d[$i].sW = $(d[i].sW)") 
end

for i=1:5 
	println("d[$i].pW = $(d[i].pW)") 
end

for i=2:5 
	println("d[$i].pw = $(flattened(d[i].pW))") 
end

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
cmaes = runEA(sys, f; center_init = ones(n), σ_init = 1.0, verbose = ReturnLevel(), monitoring = true);
ipop = runEA(sys, rsys, f; center_init = ones(n), σ_init = 1.0, verbose = RestartLevel(), monitoring = true);

# for levy and rosenbrock
cmaes = runEA(sys, f; center_init = zeros(n), σ_init = 1.0, verbose = ReturnLevel());
ipop = runEA(sys, rsys, f; center_init = zeros(n), σ_init = 1.0, verbose = RestartLevel());

# Write output to file
outputPath = "/Users/mwineberg/Documents/2 Work/Projects/CMAES_DE"

rInfo = cmaes.runInfo;
write_selectionsource(sys, rInfo; path = outputPath)

rInfo = ipop.runInfo;
write_selectionsource(sys, rsys, rInfo; path = outputPath)

rand_shell(n::Int; radius = sqrt(n)) = rand_shell(zeros(n); radius = radius)

# radius default is the same distance as euclidian([0,0,..,0], [1,1,..,1])
function rand_shell(center::Vector; radius = sqrt(length(center)))
	ci = rand(length(center))
	sum_ci_sq = sum(map((x)->x^2, ci))
	ci = radius / sqrt(sum_ci_sq) * ci + center
end

# Example of running mulitple experiments and writing ReturnInfo
function runexpr(exprName::String; reps = 20, outputPath = "")
	prefixNames = ["fn", "dim", "elitism", "center", "run"]
	firstTime = true
	for n = [25, 50]
		testFn = generatetests(n, 0.0; ε = 1.0e-5)
		fn_name = [:rastrigin, :levy, :elliptical, :ackley, :griewank]
		for name in fn_name
			f = testFn[name]
			r_UnitShell = rand_shell(optimum(f))
			for includeCenter = [false,true], elitism = [false,true]
				sys = CMAES_System(n, f; maxEvals = 1_000_000, includeCenter = includeCenter, elitism = elitism)
				rsys = CMAES_Restart(; η = 2)
				deg = 5.0
				for expr = 1:reps
					prefixValues = [name, n, elitism, includeCenter, expr]
					println("\n\n-------------------------------------------------------------------------------------------------------------")
					println("Fn = $name, dim = $n, elitism = $elitism, includeCenter = $includeCenter, run = $(expr)/$reps")
					ipop = runEA(sys, rsys, f; center_init = r_UnitShell, σ_init = 1.0, verbose = RestartLevel(), monitoring = false)
					write_final(ipop; prefixNames = prefixNames, prefixValues = prefixValues, 
						              initialize = firstTime, path = outputPath, fileName = "ipop+_final$exprName")
					# write_run(ipop, sys, f; prefixNames = prefixNames, prefixValues = prefixValues,  
					#    			         initialize = firstTime, path = outputPath, fileName = "ipop+_run$exprName")
					firstTime = false
				end
			end
		end
	end
end

expr_path = "/Users/mwineberg/Dropbox/2 Work/Students/Graduate Students (Advisor)/Opawale, Samuel/Programs/Experiments/ipop"
runexpr("<single run all>", reps = 1, outputPath = expr_path)

# Look inside final state of the first rep
rState = ipop.state;
r1State = first(rState);
r1popn = population(r1State, :post);

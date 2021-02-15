rand_shell(n::Int; radius = sqrt(n)) = rand_shell(zeros(n); radius = radius)

# radius default is the same distance as euclidian([0,0,..,0], [1,1,..,1])
function rand_shell(center::Vector; radius = sqrt(length(center)))
	ci = rand(length(center))
	sum_ci_sq = sum(map((x)->x^2, ci))
	ci = radius / sqrt(sum_ci_sq) * ci + center
end

# Example of running mulitple experiments and writing ReturnInfo
# IMPORTANT: never turn monitored to false. The dualcenter algorithm now uses the monitor data in it's model
function runexpr(exprName::String; reps = 20, outputPath = "", summary = true, monitored = true)
	#prefixNames = ["fn", "dim", "elitism", "ctr", "run"]
	prefixNames = ["fn", "dim","run"]
	firstTime = true
	for n = [40]
		testFn = generatetests(n, 0.0; ε = 10e-8)
		fn_name  = [:rastrigin]
		for name in fn_name
			(name == :rastrigin) ? max_fe = 3 * 10000 * n : max_fe = 10000n
			f = testFn[name]
			r_UnitShell = rand_shell(optimum(f))
			for includeCenter = [false], elitism = [false]
				sys = CMAES_System(n, f; maxEvals = max_fe, includeCenter = includeCenter, elitism = elitism)
				rsys = CMAES_Restart(; η = 2.0)
				deg = 5.0
				for expr = 1:reps
					#prefixValues = [name, n, elitism, includeCenter, expr]
					prefixValues = [name, n, expr]
					println("\n\n-------------------------------------------------------------------------------------------------------------")
					println("Fn = ",name," dim = ",n," elitism = ",elitism," includeCenter = ",includeCenter, " run = ",expr,"/",reps)
					ipop = runEA(sys, rsys, f; center_init = r_UnitShell, σ_init = 1.0, verbose = RestartLevel(), monitoring = monitored)

					if summary
						write_final(ipop; prefixNames = prefixNames, prefixValues = prefixValues,
						            	  initialize = firstTime, path = outputPath, fileName = "Rastrigin-40D-(3)10000n-rs1505-test$exprName")
					end

					#=if monitored
						write_run(ipop, sys, f; prefixNames = prefixNames, prefixValues = prefixValues,
					   			    	     	initialize = firstTime, path = outputPath, fileName = "allfns_euw_dual0220_run$exprName", sep = ",")
					end=#

					firstTime = false
				end
			end
		end
	end
end

expr_path = "$(base_path)/Experiments/fixedbudgetruns"
runexpr("#test", reps = 5, outputPath = expr_path, monitored = true)

function de_test(fit; popn_size = 100, max_evals = 100_000, verbose = ReturnLevel(), selection = slotselection, segCount = 1,
                   monitoring = false, restart = false, template = [:classic, :classic, :classic_direct, :direct], 
                   donorSel = :uniform, parentSel = :uniform, localDonor = 1.0, CR = 0.5)
  de_sys = DE_System(fit, max_evals; popnSize = popn_size, selection = selection, segCount = segCount, 
  	                 donorSel = donorSel, parentSel = parentSel, localDonor = localDonor, CR = CR)
  if restart
  	r_sys = DE_Restart(; η = 2)
  	runEA(de_sys, r_sys, fit; verbose = verbose, monitoring = monitoring);
  else
  	runEA(de_sys, fit; verbose = verbose, monitoring = monitoring);
  end
end

rand_shell(n::Int; radius = sqrt(n)) = rand_shell(zeros(n); radius = radius)

# radius default is the same distance as euclidian([0,0,..,0], [1,1,..,1])
function rand_shell(center::Vector; radius = sqrt(length(center)))
	ci = rand(length(center))
	sum_ci_sq = sum(map((x)->x^2, ci))
	ci = radius / sqrt(sum_ci_sq) * ci + center
end

# Example of running mulitple experiments and writing ReturnInfo
function runexpr(exprName::String; reps = 20, 
	             outputPath = "/Users/mwineberg/Dropbox/2 Work/Students/Graduate Students (Advisor)/Opawale, Samuel/Programs/Experiments/de")
	prefixNames = ["fn", "dim", "deg", "template", "sel", "restart", "run"]
	firstTime = true
	fn_name = [:rastrigin, :levy, :elliptical, :griewank]
	for n = [10], deg = [0.0, 5.0]
		testFn = boundedtests(n, deg; ε = 1.0e-5)
		for name in fn_name
			f = testFn[name]
			r_UnitShell = rand_shell(optimum(f))
			for template = [:best, :classic], sel = [slotselection, truncationselection], restart = [true], expr = 1:reps
				prefixValues = [name, n, deg, template, sel, restart, expr]
				println("\n\n-------------------------------------------------------------------------------------------------------------")
				println("Fn = $name, dim = $n, deg = $deg, template = $template, selection = $sel, restart = $restart, run = $(expr)/$reps")
				de = de_test(f, max_evals = 100_000, popn_size = 32, template = [template], 
								selection = sel);
				# de = de_test(f, max_evals = 500_000, popn_size = :default, restart = restart, template = [template], 
				# 				selection = sel, verbose = RestartLevel());
				write_final(de; prefixNames = prefixNames, prefixValues = prefixValues, 
					              initialize = firstTime, path = outputPath, fileName = "DE_final$exprName")
				# write_run(ipop, sys, f; prefixNames = prefixNames, prefixValues = prefixValues,  
				#    			         initialize = firstTime, path = outputPath, fileName = "ipop+_run$exprName")
				firstTime = false
			end
		end
	end
end

# runexpr("<n=5_exp_w_norestart_2>", reps = 20)
expr_path = "$(base_path)/Experiments"
runexpr("<n=10_exp_w_norestart_2>", reps = 1, outputPath = expr_path)


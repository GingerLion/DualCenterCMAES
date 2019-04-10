# rosenbrock() = Real
# rosenbrock(x::Vector{Real}) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# rosenbrock(x::Vector{Float64}) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# fit = DE_Fitness(rosenbrock, 2, :min, 0.0, 1e-10, sqrt(3), -sqrt(3))


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

btestFn3 = boundedtests(3)
btestFn3 = boundedtests(3, 5.0)		# d = 3, rotated by 5 degress
f = btestFn3[:rastrigin]
f = btestFn3[:schwefel]
f = btestFn3[:sphere]
f = btestFn3[:ackley]
f = btestFn3[:levy]
f = btestFn3[:griewank]

btestFn20 = boundedtests(20)   #d=20
btestFn20 = boundedtests(20, 5.0)
f = btestFn20[:griewank]
f = btestFn20[:rastrigin]
f = btestFn20[:sphere]
f = btestFn20[:elliptical]

btestFn10 = boundedtests(10; ε = 1.0e-5)   #d=20
btestFn10 = boundedtests(10, 5.0; ε = 1.0e-5)
f = btestFn10[:griewank]
f = btestFn10[:rastrigin]
f = btestFn10[:sphere]


btestFn5 = boundedtests(5; ε = 1.0e-5)   #d=20
btestFn5 = boundedtests(5, 5.0; ε = 1.0e-5)
f = btestFn5[:rastrigin]
f = btestFn5[:griewank]

btestFn2 = boundedtests(2)
btestFn2 = boundedtests(2, 5.0)		# d = 2, rotated by 5 degress
f = btestFn2[:rastrigin]
f = btestFn2[:schwefel]
f = btestFn2[:sphere]
f = btestFn2[:ackley]
f = btestFn2[:levy]
f = btestFn2[:griewank]

btestFn5 = boundedtests(5)
f = btestFn5[:levy]
f = btestFn5[:rastrigin]

btestFn10 = boundedtests(10)
f = btestFn10[:rastrigin]
f = btestFn10[:levy]

de = de_test(f, max_evals = 20000, verbose = true);

de = de_test(f, max_evals = 20000, template =[:best]);
de = de_test(f, max_evals = 20000, template =[:direct]);
de = de_test(f, max_evals = 20000, template =[:classic; :classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, template =[:best]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 5, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 50, template =[:classic]);

de = de_test(f, max_evals = 20000, template =[:classic; :classic]);

de = de_test(f, max_evals = 20000, template =[:classic], donorSel = :rank, parentSel = :rank);
de = de_test(f, max_evals = 20000, template =[:classic], selection = truncationselection);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 2, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 5, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 10, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 20, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 25, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 50, template =[:classic]);

de = de_test(f, max_evals = 20000, template =[:classic], donorSel = :rank);
de = de_test(f, max_evals = 20000, template =[:classic], selection = truncationselection);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 2, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 5, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 10, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 20, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 25, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 50, template =[:classic]);

de = de_test(f, max_evals = 20000, template =[:classic]);
de = de_test(f, max_evals = 20000, template =[:best]);
de = de_test(f, max_evals = 20000, template =[:classic], selection = truncationselection);
de = de_test(f, max_evals = 20000, template =[:best], selection = truncationselection);
de = de_test(f, max_evals = 20000, template =[:best]);
de = de_test(f, max_evals = 20000, template =[:best], parentSel = :rank);
de = de_test(f, max_evals = 20000, template =[:best], selection = truncationselection);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 2, template =[:best]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 5, template =[:best]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 10, template =[:best]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 20, template =[:best]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 25, template =[:best]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 50, template =[:best]);


de = de_test(f, max_evals = 20000, popn_size = 30, template =[:classic]);
de = de_test(f, max_evals = 20000, selection = truncationselection, segCount = 5, template =[:classic]);


de = de_test(f, max_evals = 20000, popn_size = 100, CR = 0.5, template =[:best], selection = truncationselection);
de = de_test(f, max_evals = 20000, template =[:classic], selection = truncationselection);

de = de_test(f, max_evals = 100000, popn_size = 256, template =[:best], selection = truncationselection);
de = de_test(f, max_evals = 100000, popn_size = :default, restart = true, template =[:best], selection = truncationselection, verbose = RestartLevel());
de = de_test(f, max_evals = 100000, popn_size = :default, restart = true, template =[:classic], localDonor = 0.25, selection = truncationselection, verbose = RestartLevel());
de = de_test(f, max_evals = 100000, popn_size = :default, restart = true, template =[:classic], localDonor = 1.0, selection = truncationselection, verbose = RestartLevel());

de = de_test(f, max_evals = 100000, popn_size = 256, template =[:best], selection = truncationselection, verbose = RunLevel());
de = de_test(f, max_evals = 100000, popn_size = :default, restart = true, template =[:best], selection = truncationselection, verbose = RunLevel());

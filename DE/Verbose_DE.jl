#----------------------------------------------------------
# State Printing
#----------------------------------------------------------
# gen::Int64
# evalCount::Int64
# contenders::Contenders
# popn::Population
# status::Symbol
# best::Tuple
# sys::DE_System
#----------------------------------------------------------

function println(state::DE_State) 
  (bestChr, bestFit) = best(state)
  println("gen = $(currentgen(state)),  local-eval = $(evals(state))")
  println("          best: fit[$bestFit]")
  println("          chr$bestChr")
#  println("popn = $(population(state))")
end


function println(state::DE_State, restart::RestartState) 
  (bestChr, bestFit) = best(state)
  stagnation = restart.stagnation
  println("\nRestart: rep = $(rep(restart)) at gen = $(currentgen(state))")
  println("         local-eval = $(evals(state)), total-eval = $(evals(restart))")
  println("          best: fit[$bestFit]")
  println("          chr$bestChr")
  conditionPrinted = false
#  msg = ["hist", "tolx", "zeroAxis", "negEigVal", "complexEigVal", "zeroEigVal", "zeroCoord", "cond"]
  msg = ["hist"]
  if any(stagnation)
    print("<")

    for i = 1:length(stagnation)
      if stagnation[i] 
        print(conditionPrinted ? ", " : "") 
        print(msg[i]) 
        conditionPrinted = true 
      end
    end

    println(">")
  end
end

#----------------------
# Printing states

function println(state::CMAES_State) 
  print("$(status(state)): ")
  print("[gen = $(gen(state))]")
  print("[σ=$(sigma(state))]")
  print("[h_σ = $(state.nModel.h_σ)]")
  print("  best = ")
  print("fit[$(state.sOffspring[:fit, 1])]")
  println(" chr$(state.sOffspring[:chr, 1])")
end

function println(state::CMAES_State, restart::RestartState) 
  stagnation = restart.stagnation
  println("\nRestart: rep = $(rep(restart)) at gen = $(currentgen(state))")
  conditionPrinted = false
  msg = ["hist", "tolx", "zeroAxis", "negEigVal", "complexEigVal", "zeroEigVal", "zeroCoord", "cond"]
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

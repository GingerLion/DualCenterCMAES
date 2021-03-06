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
  stagnation_ = restart.stagnation_
  println("\nRestart: rep = ",rep(restart)," at gen = ",currentgen(state))
  println(firstRequest(state)," system requested to restart first.")
  if !isempty(restart.stagReason)
    conditionPrinted = false
    print("<")
    for i in restart.stagReason
        conditionPrinted ? print(",") : print("")
        print(i)
        conditionPrinted = true
    end
    print("> : normal system\n")
  end
  if !isempty(restart.stagReason_)
      conditionPrinted = false
      print("<")
      for i in restart.stagReason_
          conditionPrinted ? print(",") : print("")
          print(i)
          conditionPrinted = true
      end
      print("> : dualcenter system\n")
  end
  #=conditionPrinted = false
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

    println("> : normal system")
end
    conditionPrinted = false
    msg = ["hist","tolx", "zeroAxis", "negEigVal", "complexEigVal", "zeroEigVal", "zeroCoord", "cond"]
  if any(stagnation_)
    print("<")

    for i = 1:length(stagnation_)
      if stagnation_[i]
        print(conditionPrinted ? ", " : "")
        print(msg[i])
        conditionPrinted = true
      end
    end

    println("> : dualcenter system")
end=#
end

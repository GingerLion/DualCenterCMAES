function DE_System(fit::DE_Fitness, maxEvals::Integer; popnSize = :default, 
                 CR = 0.5, diffWeight = 0.85, initType = :unif,
                 selection = slotselection, segCount = 1, randomizeSlots = false,
                 donorSel = :uniform, donorSlope = 2, parentSel = :uniform, parentSlope = 2, 
                 localDonor = 1.0, template = [:classic])  
  initFn = PopnInitializer(fit)[initType]
  chrLength = fit.dimension

  popnSize = ((popnSize == :default) ? 4 + floor(Int, 3 * log(chrLength)) : popnSize)

  sp = DE_SamplingParms(CR, diffWeight, donorSel, donorSlope, parentSel, parentSlope, localDonor, template)

  samplers = samplertemplate(SlotSamplers(popnSize, sp), template)
  evalsPerGen = popnSize * length(samplers)

  DE_System(maxEvals, popnSize, chrLength, selection, segCount, randomizeSlots, initFn, sp, samplers, evalsPerGen)
end

function DE_System(de::DE_System, fit::DE_Fitness, restart::RestartState, verbose::Verbose)
  popnSize = nextpopnsize(de, restart)
  samplers = samplertemplate(SlotSamplers(popnSize, de.rParms), de.rParms.template)
  evalsPerGen = popnSize * length(samplers)

  if atlevel(verbose, RestartLevel)
    println("new popnSize = $popnSize")
  end

 DE_System(de.maxEvals, popnSize, de.chrLength, de.selection, de.segCount, 
           de.randomizeSlots, de.initialize, de.rParms, samplers, evalsPerGen)
end

evalsPerGen(sys::DE_System) = sys.evalsPerGen
popnsize(sys::DE_System) = sys.popnSize
chrlength(sys::DE_System) = sys.chrLength

function runEA(sys::DE_System, fit::DE_Fitness; monitoring = false, verbose = NotVerbose())
  runInfo = newmonitor(monitoring)
  state = DE_State(sys, fit, runInfo, verbose)
  returnInfo = ReturnInfo(state, runInfo)
  runEA(state, fit, runInfo, returnInfo, verbose)
end

function runEA(sys::DE_System, rsys::DE_RestartBase, fit::DE_Fitness; monitoring = false, verbose = NotVerbose())
  runInfo = newmonitor(monitoring)
  state = DE_State(sys, fit, runInfo, verbose)
  returnInfo = ReReturnInfo(state, runInfo)
  rsys = DE_RestartFull(rsys, state)
  restart = RestartState(rsys, state)
  runEA(state, restart, fit, runInfo, returnInfo, verbose)
end

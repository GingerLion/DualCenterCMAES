#-----------------------------------------------------------------------
#  CMAES_System <: System   (default implementation as CMA_ES_System <: CMAES_System)
#-----------------------------------------------------------------------
#   maxEvals::Int
#   rParms::Reproduction_Parms  # constants used during  reproduction (aka Model_Parms)
#   sParms::Selection_Parms     # constants used during selection
#-----------------------------------------------------------------------

# -------------------------------
# External Constructor

function CMAES_System(n::Integer, f::RealFitness;
                      maxEvals = 10_000, modelParms = :default,
                      segmentCount = 1, beginBinding = :mostFit, loaded = :mostFit,
                      includeCenter = false, elitism = false, λ = 4 + floor(Int, 3 * log(n)), η = 1)
  modelParms = (modelParms == :default) ? Model_Parms(n, f; λ = λ) : modelParms
  selectionParms = Selection_Parms(n, modelParms.μ, modelParms.λ, f.direction,
                                   loaded, beginBinding, segmentCount, includeCenter, elitism, η)
  CMAES_System(maxEvals, modelParms, selectionParms)
end

function CMAES_System(sys::CMAES_System, f::RealFitness, restart::RestartState, verbose::Verbose)
  (μ, λ) = nextpopnsize(sys, restart)

  if atlevel(verbose, RestartLevel)
    println("new μ = $μ, new λ = $λ")
  end

  n = sys.rParms.N
  sel = sys.sParms
  rParm = Model_Parms(n, f; μ = μ, λ = λ)
  sel = Selection_Parms(n, μ, λ, sel.direction, sel.loaded, sel.beginBinding,
                        sel.segmentCount, sel.includeCenter, sel.elitism, sel.η)
  CMAES_System(sys.maxEvals, rParm, sel)
end


# -------------------------------
# General Functions

selection_parms(sys::CMAES_System) = sys.sParms

includecenter(sys::CMAES_System) = sys.sParms.includeCenter
elitism(sys::CMAES_System) = sys.sParms.elitism
segmentcount(sys::CMAES_System) = sys.sParms.segmentCount

η(sys::CMAES_System) = sys.sParms.η
mu(sys::CMAES_System) = sys.sParms.μ
lambda(sys::CMAES_System) = sys.sParms.λ
weights(sys::CMAES_System) = sys.rParms.w
popnsize(sys::CMAES_System) = (mu(sys), lambda(sys))
evalsPerGen(sys::CMAES_System) = lambda(sys) + (includecenter(sys) ? 1 : 0)
evalsPerGen_(sys::CMAES_System) = 2 * lambda(sys) + (includecenter(sys) ? 1 : 0)

function runEA(sys::CMAES_System, f::RealFitness; center_init = ones(n), σ_init = 1.0,
                   monitoring = false, verbose = NotVerbose())
  runInfo = newmonitor(monitoring, sys)
  state = CMAES_State(sys, f, center_init, σ_init, runInfo, verbose)
  returnInfo = ReturnInfo(state,runInfo)
  runEA(state, f, runInfo, returnInfo, verbose)
end

function runEA(sys::CMAES_System, rsys::CMAES_RestartBase, f::RealFitness; center_init = ones(n), σ_init = 1.0,
                   monitoring = false, verbose = NotVerbose())
  runInfo = newmonitor(monitoring, sys)
  state = CMAES_State(sys, f, center_init, σ_init, runInfo, verbose)
  returnInfo = ReReturnInfo(state, runInfo)
  rsys = CMAES_RestartFull(rsys, state, center_init, σ_init)
  restart = RestartState(rsys, state)
  runEA(state, restart, f, runInfo, returnInfo, verbose)
end

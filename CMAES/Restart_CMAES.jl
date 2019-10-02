
#-----------------------------------------------------------------------
#  Restart System
#-----------------------------------------------------------------------
# struct CMAES_Restart <: Restart
#   testCount::Int
#   η::Int,
#   η_μ::Int,
#   η_λ::Int,
#   initλ::Int,
#   tol_f::Float64          # const   CMAES only (probably)
#   tol_x::Float64          # const   CMAES only
#-----------------------------------------------------------------------
# mutable struct RestartState
#   rep::Int
#   totalEvals::Integer
#   bfHist::BestFitHistory
#   stagnation::Vector{Bool}
#   parms::Restart
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Outer Constructors
#-----------------------------------------------------------------------

function CMAES_Restart(; historyWindow = :default, ignoreStagnation = false,
	                   η = 2.0, η_μ = η, η_λ = η, tol_f = 10.0^-12)
  CMAES_RestartBase(ignoreStagnation, historyWindow, η, η_μ, η_λ, tol_f)
end

function CMAES_RestartFull(rsys::CMAES_RestartBase, state::CMAES_State, center_init::Vector,  σ_init::Float64)
  tol_x = rsys.tol_f * sigma(state)
  tol_x_ = rsys.tol_f * sigma_(state)
  initLambda = lambda(state)
  #σ_mult = 2*σ_init
  hw = ((rsys.historyWindow == :default) ? 10 + ceil(Int, 30 * dimensions(state) / lambda(state))
                                         : rsys.historyWindow)
  CMAES_RestartFull(rsys.ignoreStagnation, hw, rsys.η, rsys.η_μ, rsys.η_λ, initLambda, rsys.tol_f, tol_x, tol_x_, center_init, σ_init)
end


#-----------------------------------------------------------------------
# General Functions - CMAES_Restart only
#-----------------------------------------------------------------------
tol_x(rsys::CMAES_Restart) = rsys.tol_x
tol_x_(rsys::CMAES_Restart) = rsys.tol_x_


#-----------------------------------------------------------------------
# General Functions - RestartState only
#-----------------------------------------------------------------------

scalemu(restart::RestartState) = restart.parms.η_μ
scalelambda(restart::RestartState) = restart.parms.η_λ
initλ(restart::RestartState) = restart.parms.initλ
initcenter(restart::RestartState) = restart.parms.center_init
initσ(restart::RestartState) = restart.parms.σ_init
tol_f(restart::RestartState) = restart.parms.tol_f
tol_x(restart::RestartState) = restart.parms.tol_x
tol_x_(restart::RestartState) = restart.parms.tol_x


#-----------------------------------------------------------------------
# RestartState & CMAES_State interactions
#-----------------------------------------------------------------------

function stagnationcriteria(state::CMAES_State, restart::RestartState)
  (model, gen) = (currentmodel(state), currentgen(state))
  [equalfunvalhist(state, restart),
   #toomuchfluctuation(state, restart),
   tolx(model, tol_x(restart)),
   noeffectaxis(model, gen),
   negeigerr(state),
   complexeigerr(state),
   zeroeigerr(state),
   noeffectcoord(model),
   conditioncov(model)]
end

function stagnationcriteria_(state::CMAES_State, restart::RestartState)
  (model_shadow, gen_shadow) = (currentmodel_(state), currentgen_(state))
  [equalfunvalhist_(state, restart),
   #toomuchfluctuation_(state, restart),
   tolx(model_shadow, tol_x_(restart)),
   noeffectaxis(model_shadow, gen_shadow),
   negeigerr_(state),
   complexeigerr_(state),
   zeroeigerr_(state),
   noeffectcoord(model_shadow),
   conditioncov(model_shadow)]
end

function nextpopnsize(sys::CMAES_System, restart::RestartState)
  (convert(Int,floor(mu(sys) * scalemu(restart))), convert(Int,floor(lambda(sys) * scalelambda(restart))))
end

function next!(restart::RestartState, state::CMAES_State)
  restart.rep += 1
  histWindow = 10 * restart.rep + ceil(Int, 30 * dimensions(state) / initλ(restart))
  restart.bfHist = BestFitHistory(histWindow, direction(state))
  restart.bfHist_ = BestFitHistory(histWindow, direction(state))
  restart.bcHist_ = BestChrHistory(10, direction(state))
  restart.bcHist_fitnesses_ = BestFitHistory(10, direction(state))
  restart.sortOrder = Vector{Int}(undef,10) # used when using sortedelitewindow
  restart.shouldRestart = false
end

#---------------------------------------------------------------------------------------------------------------------------------
# mutable struct CMAES_State <: State  (default implementation as CMA_ES_State <: CMAES_State)
#---------------------------------------------------------------------------------------------------------------------------------
#   gen::Integer
#	gen_shadow::Integer		    #shadow shadow
#   evalCount::Integer
#	evalCount_shadow::Integer		#shadow shadow
#   sys::CMAES_System           # constants used during the running of the system (includes reproduction and selection parameters)
#   pModel::CMAES_Model         # model from previous gen, used at the beginning of reproduction step
#   pModel_shadow::CMAES_Model		# shadowged model from previous gen
#   pOffspring::Population      # popnsize = μ    offspring from previous gen  (used for elitism)
#	pOffspring_shadow::Population		# popnsize = μ	  shadowged offspring from previous gen (used for elitism)
#   pE::SphericalNoise          # popnsize = μ    E from previous gen         (currenlty used for elitism) shadowGED doesnt matter nE is all same
#   pW::ShapedNoise             # popnsize = μ    W from previous gen         (used for elitism)
#	pW_shadow::ShapedNoise         # popnsize = μ    shadowged W from previous gen         (used for elitism)
#   nE::SphericalNoise          # popnsize = λ    generated noise used to create samples
#   nW::ShapedNoise             # popnsize = λ    shaped noise used to create samples
#   nW_shadow::ShapedNoise         # popnsize = λ    shadowged shaped noise used to create samples
#   nOffspring::Population      # popnsize = λ    samples used to create new model after selection
#	nOffspring_shadow::Population  # popnsize = λ    shadowged samples used to create new model after selection
#   sOffspring::Population      # popnsize = μ;   sorted and truncated (selection) changing sample distribution for model updating
#	sOffspring_shadow::Population      # popnsize = μ;   shadowged sorted and truncated (selection) changing sample distribution for model updating
#   sW::ShapedNoise             # popnsize = μ;   sorted and truncated (selection) changing sample distribution for model updating
#   sW_shadow::ShapedNoise             # popnsize = μ;   shadowged sorted and truncated (selection) changing sample distribution for model updating
#   nModel::CMAES_Model         # model after the reproductive step complete
#	nModel_shadow::CMAES_Model     # shadowged model after the reproductive step complete
#   status::Symbol
#	status_shadow::Symbol			# shadowged status
#   best::Tuple
#	best_shadow::Tuple				# shadowged best chromosome
#	firstRequest::Symbol		# records the first system that wishes to restart
#---------------------------------------------------------------------------------------------------------------------------------

function setup!(state::CMAES_State, model::CMAES_Model, sys::CMAES_System, f::RealFitness)
	(n, λ, μ, direction) = (sys.rParms.N, sys.sParms.λ, sys.sParms.μ, sys.sParms.direction)
	state.nModel = model
	state.nModel_shadow = deepcopy(model)
	state.sys = sys
	state.gen = 0
	state.gen_shadow = state.gen
	state.evalCount = 0
	state.evalCount_shadow = 0
	state.nE = ZeroNoise(SphericalNoise, n, λ)
	state.nW = ZeroNoise(ShapedNoise, n, λ)
	state.nW_shadow = deepcopy(state.nW)
	state.sW = ZeroNoise(ShapedNoise, n, μ)
	state.sW_shadow = deepcopy(state.sW)
	state.pW = deepcopy(state.nW)
	state.pW_shadow = deepcopy(state.pW)
	state.sOffspring = SortedPopulation(RegularPopulation(center(model), μ, f.objFn; direction = direction))
	state.sOffspring_shadow = deepcopy(state.sOffspring)
	state.pOffspring = deepcopy(state.sOffspring)
	state.pOffspring_shadow = deepcopy(state.pOffspring)
	state.best = best(state.sOffspring)
	state.best_shadow = best(state.sOffspring_shadow)
	state.firstRequest = :none
	evolvable!(state)
	evolvable!_(state)
end


# -------------------------------
# External Constructors

function CMAES_State(sys::CMAES_System, f::RealFitness,
	                 center_init::Vector, σ_init::Float64,
	                 runInfo::Monitor, verbose::Verbose = NotVerbose())
  model = CMAES_Model(sys.rParms, copy(center_init), σ_init)
  CMAES_State(model, sys, f, runInfo, verbose)
end

function CMAES_State(sys::CMAES_System, f::RealFitness, restart::RestartState,
                   runInfo::Monitor, verbose::Verbose = NotVerbose())
  model = CMAES_Model(sys.rParms, copy(initcenter(restart)), initσ(restart))
  CMAES_State(model, sys, f, restart, runInfo, verbose)
end



# -------------------------------
# General functions

parms(c::CMAES_State) = (c.sys, c.sys.sParms, c.sys.rParms)
sourceparms(c::CMAES_State) = c.sourceParms
errorstate(c::CMAES_State) = (c.status == :error)
error!(c::CMAES_State, errorCode::Symbol = :error) = (c.status = errorCode)
error!_(c::CMAES_State, errorCode::Symbol = :error) = (c.status_shadow = errorCode)
currentmodel(c::CMAES_State) = c.nModel
mu(c::CMAES_State) = c.nModel.parms.μ
lambda(c::CMAES_State) = c.nModel.parms.λ
sigma(c::CMAES_State) = c.nModel.σ
center(c::CMAES_State) = c.nModel.center
covar(c::CMAES_State) = c.nModel.C
parms(sys::CMAES_System) = (sys.sParms, sys.rParms)
segmentcount(c::CMAES_State) = c.sys.sParms.segmentCount
initializing(c::CMAES_State) = (c.gen == 0)
pW(c::CMAES_State) = c.pW
sW(c::CMAES_State) = c.sW
firstRequest(c::CMAES_State) = c.firstRequest

#General shadowed Functions
currentmodel_(c::CMAES_State) = c.nModel_shadow
#mu_(c::CMAES_State) = c.nModel_shadow.parms.μ
#lambda_(c::CMAES_State) = c.nModel_shadow.parms.λ
sigma_(c::CMAES_State) = c.nModel_shadow.σ
center_(c::CMAES_State) = c.nModel_shadow.center
covar_(c::CMAES_State) = c.nModel_shadow.C
initializing_(c::CMAES_State) = (c.gen_shadow == 0)
pW_(c::CMAES_State) = c.pW_shadow
sW_(c::CMAES_State) = c.sW_shadow
population_shadow(c::CMAES_State) = c.sOffspring_shadow



# note: currently during evolution parents are in :post while :parents holds parents of the parents
#       after evolution parents are in :parents and :post holds the post selection population
#       (to be changed so that both during and after evolution :parents holds parents
#        while for :post: before truncation it will hold the parents;
#                         after truncation it will hold the truncated offspring)
function population(state::CMAES_State, selection::Symbol = :post)
  if selection == :pre
  	state.nOffspring
  elseif selection == :post
  	state.sOffspring
  elseif selection == :parents
    state.pOffspring
  else
  	error("selection choice should be either :pre, :post or :parents; instead it is $selection")
  end
end
#for shadowged populations
function population_(state::CMAES_State, selection::Symbol = :post)
  if selection == :pre
  	state.nOffspring_shadow
  elseif selection == :post
  	state.sOffspring_shadow
  elseif selection == :parents
    state.pOffspring_shadow
  else
  	error("selection choice should be either :pre, :post or :parents; instead it is $selection")
  end
end

function model(state::CMAES_State, selection::Symbol = :post)
  if selection == :pre
    state.pModel
  elseif selection == :post
    state.nModel
  else
    error("selection choice should be either :pre or :post; instead it is $selection")
  end
end
#for shadowged models
function model_(state::CMAES_State, selection::Symbol = :post)
  if selection == :pre
    state.pModel_shadow
  elseif selection == :post
    state.nModel_shadow
  else
    error("selection choice should be either :pre or :post; instead it is $selection")
  end
end


sortorder(c::CMAES_State) = SortStructure(population_(c, :post)) # changed to population_ to view the shadow's selection source
sortstructure(c::CMAES_State) = SortStructure(population_(c, :post))   # synonym to sortorder()

negeigerr(state::CMAES_State) = (status(state) == :negativeEigenError)
complexeigerr(state::CMAES_State) = (status(state) == :complexEigenError)
zeroeigerr(state::CMAES_State) = (status(state) == :zeroEigenError)
eigenerror(state::CMAES_State) = (negeigerr(state) || complexeigerr(state) || zeroeigerr(state))
negeigerr_(state::CMAES_State) = (status_(state) == :negativeEigenError)
complexeigerr_(state::CMAES_State) = (status_(state) == :complexEigenError)
zeroeigerr_(state::CMAES_State) = (status_(state) == :zeroEigenError)
eigenerror_(state::CMAES_State) = (negeigerr_(state) || complexeigerr_(state) || zeroeigerr_(state))

function invsqrtC!(c::CMAES_State, m::CMAES_Model; shadow = false)
  if !shadow
	  try
	  	invsqrtC!(m)
	  catch
	    error!(c, :zeroEigenError)
	  end
  else
	  try
		  invsqrtC!(m)
	  catch
		  error!_(c, :zeroEigenError)
	  end
  end
end

function eigendecomp!(c::CMAES_State, m::CMAES_Model; shadow = false)
  if !shadow
	  try
	  	eigendecomp!(m)
	  catch
	    error!(c, :complexEigenError)
	  end
  else
	  try
		  eigendecomp!(m)
	  catch
		  error!_(c, :complexEigenError)
	  end
  end
end

function sqrtC!(c::CMAES_State, m::CMAES_Model; shadow = false)
	if !shadow
	  try
	  	sqrtC!(m)
	  catch
	    error!(c, :negativeEigenError)
	  end
  	else
		try
			sqrtC!(m)
		catch
			error!_(c, :negativeEigenError)
		end
	end
end


# -------------------------------
# Restart related functions

evolvable(state::CMAES_State, restart::RestartState) = (evolvable(state) || eigenerror(state))
evolvable_(state::CMAES_State, restart::RestartState) = (evolvable_(state) || eigenerror_(state))

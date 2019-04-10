#----------------------
# mutable struct CMAES_RunInfo <: RunInfo
#   data::Dict
#   uq::Tuple{Int64,Float64}
#   med::Tuple{Int64,Float64}
#   lq::Tuple{Int64,Float64}
# 	sourceValues::SelectionSourceParms
# 	allSourceValues::SelectionSourceParms
#----------------------
#   CMAES_RunInfo(w::Weights)
#----------------------

newmonitor(monitoring::Bool, sys::CMAES_System) =  (monitoring ? CMAES_RunInfo(sys) : NoMonitor())

sourcevalues(runInfo::RunInfo) = runInfo.sourceValues
allsourcevalues(runInfo::RunInfo) = runInfo.allSourceValues

function update!(nm::NoMonitor, sys::CMAES_System) end

function update!(runInfo::RunInfo, sys::CMAES_System)
	w = weights(sys)
    runInfo.uq = locquantile(0.75, w)
    runInfo.med = locquantile(0.5, w)
    runInfo.lq = locquantile(0.25, w)
    runInfo.sourceValues = SelectionSourceParms(sys)
    runInfo.allSourceValues = AllSourceParms(sys)      # for all-mirror
end

#----------------------
# RunInfo - Updates

function monitor!(nm::NoMonitor, state::CMAES_State) end

function monitor!(runInfo::RunInfo, state::CMAES_State, f::RealFitness)
  # system models
  preModel  = model(state, :pre)
  postModel = model(state, :post)
  w = weights(preModel)

  # bugged system models
  preModel_bug = model_(state, :pre)
  postModel_bug = model_(state, :post)
  w_bug = weights(preModel_bug)

  # main system state
  parents   = population(state, :parents)
  offspring = population(state, :pre)
  selected  = population(state, :post)
  state_center = centerpopn(preModel, f)
  covMatrix = covar(postModel)

  # bugged main system state
  parents_bug = population_(state, :parents)
  offspring_bug = population_(state, :pre)
  selected_bug = population_(state, :post)
  state_center_bug = centerpopn(preModel_bug, f)
  covMatrix_bug = covar(postModel_bug)

  # general information
  runInfo[:best]    	= best(state)
  runInfo[:covar]   	= covMatrix
  runInfo[:cond]   		= cond(covMatrix)
  runInfo[:σ]       	= sigma(state)
  runInfo[:gen]			= gen(state)
  runInfo[:l_evals] 	= evals(state)

  # general information - bugged
  runInfo[:best_bug] 	= best_(state)
  runInfo[:covar_bug]   = covar_(state)
  runInfo[:cond_bug]	= cond(covMatrix_bug)
  runInfo[:σ_bug]		= sigma_(state)
  runInfo[:gen_bug]		= gen_(state)
  runInfo[:l_evals_bug]   = evals_(state)

  # system state information
  runInfo[:center]  	= (state_center[:chr, 1], state_center[:fit, 1])
  runInfo[:fitsummary]	= fitsummary(selected, runInfo, w)
  src = SelectionSource(sourcevalues(runInfo), state)
  runInfo[:source] = deepcopy(src)
  #online update of model parms
  #update!(src, postModel, system(state))
  # system state information - bugged
  runInfo[:best_fit] = bestfitness(state)
  runInfo[:best_bug_fit] = bestfitness_(state)
  runInfo[:center_bug]  = (state_center_bug[:chr, 1], state_center_bug[:fit, 1])
  #   Mirroring state
  #   note: currently does not include mirror of covar, signma or paths,
  #		    i.e. only center is mirrored from the model
  #   note: different calculation for center than state-center

  # offspring only mirror
  (offspring_selected, offs_sortOrder) = truncation(offspring, mu(state))
  offspring_center = center(offspring_selected, w)
  runInfo[:center_fit] = state_center[:fit,1]
  runInfo[:center_fit_bug] = state_center_bug[:fit,1]

  runInfo[:offspring_center] = (offspring_center, objfn(f)(offspring_center))
  runInfo[:offs_fitsummary]	 = fitsummary(offspring_selected,runInfo,  w)

  # all (parents + center + offspring) mirror
  #(all_selected, all_sortOrder) = truncation(pcat(parents, state_center, offspring), mu(state))
  #all_center = center(all_selected, w)

  #runInfo[:all_center]		 = (all_center, objfn(f)(all_center))
  #runInfo[:all_fitsummary]   = fitsummary(all_selected, runInfo, w)
  #runInfo[:all_source]       = SelectionSource(allsourcevalues(runInfo), all_sortOrder)

  #eigenvalue & eigenvector study
  if !eigenerror(state)
	  #NON-BUGGED
	  eigvals = eigendecomp(postModel)[1]
	  eigvecs = eigendecomp(postModel)[2]
	  #runInfo[:Dinv] = diag(inv(diagm(0 => eigvals)),0)
	  runInfo[:eig_vals] = eigvals
	  runInfo[:eig_vecs] = eigvecs
	  runInfo[:BD] = *(eigvecs,eigvals)
	  runInfo[:eig_min] = minimum(eigvals)
	  runInfo[:eig_max] = maximum(eigvals)
	  #BUGGED
	  eigvals_bug = eigendecomp(postModel_bug)[1]
	  eigvecs_bug = eigendecomp(postModel_bug)[2]
	  #runInfo[:Dinv_bug] = diag(inv(diagm(0 => eigvals_bug)),0)
	  runInfo[:eig_vals_bug] = eigvals_bug
	  runInfo[:eig_vecs_bug] = eigvecs_bug
	  runInfo[:BD_bug] = *(eigvecs_bug,eigvals_bug)
	  runInfo[:eig_min_bug] = minimum(eigvals_bug)
	  runInfo[:eig_max_bug] = maximum(eigvals_bug)
  else
	  #NON-BUGGED
	  runInfo[:eig_vals] = fill(NaN,chrlength(postModel))
	  runInfo[:eig_vecs] = fill(NaN,(chrlength(postModel),chrlength(postModel)))
	  #runInfo[:Dinv] = fill(NaN,chrlength(postModel))
	  runInfo[:BD] = fill(NaN,(chrlength(postModel)))
	  runInfo[:eig_min] = NaN
	  runInfo[:eig_max] = NaN
	  #BUGGED
	  runInfo[:eig_vals_bug] = fill(NaN,chrlength(postModel_bug))
	  runInfo[:eig_vecs_bug] = fill(NaN,(chrlength(postModel_bug),chrlength(postModel_bug)))
	  #runInfo[:Dinv] = fill(NaN,chrlength(postModel))
	  runInfo[:BD_bug] = fill(NaN,(chrlength(postModel_bug)))
	  runInfo[:eig_min_bug] = NaN
	  runInfo[:eig_max_bug] = NaN

  end

	#model parms
	runInfo[:μ_eff] = μ_eff(postModel)
	runInfo[:d_σ] = d_σ(postModel)
	runInfo[:chi_mean] = chi_mean(postModel)
	runInfo[:h_σ] = h_σ(postModel)
	runInfo[:c_1] = c_1(postModel)
	runInfo[:c_σ] = c_σ(postModel)
	runInfo[:c_μ] = c_μ(postModel)
	runInfo[:c_c] = c_c(postModel)
	#bugged model parms
	#not going to bother

	#paths
	runInfo[:p_σ] = p_σ(postModel)
	runInfo[:p_c] = p_c(postModel)
	runInfo[:y] = weightedavg(state.sW)
	runInfo[:invsqrtC] = invsqrtC(postModel)
	#update equations in parts
	runInfo[:p_σ_part1] = p_σ_part1(postModel)
	runInfo[:p_σ_part2] = p_σ_part2(postModel)
	runInfo[:p_σ_part3] = p_σ_part3(postModel,state.sW)

	runInfo[:σ_part1] = σ_part1(postModel)
	runInfo[:σ_part2] = σ_part2(postModel)
	runInfo[:σ_part3] = σ_part3(postModel)
	#runInfo[:σ_part4] = σ_part4(postModel)
	#runInfo[:σ_part5] = σ_part5(postModel)
	runInfo[:σ_part6] = σ_part6(postModel)

	runInfo[:p_c_part1] = p_c_part1(postModel)
	runInfo[:p_c_part2] = p_c_part2(postModel)
	runInfo[:p_c_part3] = p_c_part3(postModel,state.sW)

	runInfo[:C_part1] = C_part1(postModel)
	runInfo[:C_part2] = C_part2(postModel)
	runInfo[:C_part3] = C_part3(postModel,state.sW)

	#paths - bugged
	runInfo[:p_σ_bug] = p_σ(postModel_bug)
	runInfo[:p_c_bug] = p_c(postModel_bug)
	runInfo[:y_bug] = weightedavg(state.sW_bug)
	runInfo[:invsqrtC_bug] = invsqrtC(postModel_bug)
	#update equations in parts
	runInfo[:p_σ_part1] = p_σ_part1(postModel)
	runInfo[:p_σ_part2] = p_σ_part2(postModel)
	runInfo[:p_σ_part3] = p_σ_part3(postModel,state.sW)

	runInfo[:σ_part1] = σ_part1(postModel)
	runInfo[:σ_part2] = σ_part2(postModel)
	runInfo[:σ_part3] = σ_part3(postModel)
	#runInfo[:σ_part4] = σ_part4(postModel)
	#runInfo[:σ_part5] = σ_part5(postModel)
	runInfo[:σ_part6] = σ_part6(postModel)

	runInfo[:p_c_part1] = p_c_part1(postModel)
	runInfo[:p_c_part2] = p_c_part2(postModel)
	runInfo[:p_c_part3] = p_c_part3(postModel,state.sW)

	runInfo[:C_part1] = C_part1(postModel)
	runInfo[:C_part2] = C_part2(postModel)
	runInfo[:C_part3] = C_part3(postModel,state.sW)

	#update equations in parts - bugged
	runInfo[:p_σ_part1_bug] = p_σ_part1(postModel_bug)
	runInfo[:p_σ_part2_bug] = p_σ_part2(postModel_bug)
	#runInfo[:p_σ_part3_bug] = p_σ_part3(postModel_bug,state.sW_bug)

	runInfo[:σ_part1_bug] = σ_part1(postModel_bug)
	runInfo[:σ_part2_bug] = σ_part2(postModel_bug)
	runInfo[:σ_part3_bug] = σ_part3(postModel_bug)
	#runInfo[:σ_part4_bug] = σ_part4(postModel_bug)
	#runInfo[:σ_part5_bug] = σ_part5(postModel_bug)
	runInfo[:σ_part6_bug] = σ_part6(postModel_bug)

	runInfo[:p_c_part1_bug] = p_c_part1(postModel_bug)
	runInfo[:p_c_part2_bug] = p_c_part2(postModel_bug)
	runInfo[:p_c_part3_bug] = p_c_part3(postModel_bug,state.sW_bug)

	runInfo[:C_part1_bug] = C_part1(postModel_bug)
	runInfo[:C_part2_bug] = C_part2(postModel_bug)
	runInfo[:C_part3_bug] = C_part3(postModel_bug,state.sW_bug)
end

function monitor!(runInfo::RunInfo, state::CMAES_State, restart::RestartState, f::RealFitness)
  monitor!(runInfo, state, f)
  runInfo[:g_evals] = evals(restart)
  runInfo[:restart]	= rep(restart)
  runInfo[:lambda]  = lambda(state)
  #runInfo[:lambda_bug] = lambda_(state)
end

#-----------------------
# Helper functions

center(popn::SortedPopulation, w::Weights) = vec(mean(popn[:chr, :], w, dims=2))
centerpopn(model::CMAES_Model, f::RealFitness) = RegularPopulation(centermember(model, objfn(f)))


#-----------------------------------------
# Fitness summary calculations

function fitsummary(popn::Population, ri::CMAES_RunInfo, w::Weights)
	fit = fitness(popn)
	(maximum(fit),
	 upperquartile(fit, ri),
	 median(fit, ri),
	 mean(fit, w),
	 lowerquartile(fit, ri),
	 minimum(fit))
end

upperquartile(values::Vector{Float64}, ri::RunInfo)  = quantile(values, ri.uq)
median(values::Vector{Float64}, ri::RunInfo)         = quantile(values, ri.med)
lowerquartile(values::Vector{Float64}, ri::RunInfo)  = quantile(values, ri.lq)

function quantile(values::Vector{Float64}, loc::Tuple{Int, Float64})
	(i, percent) = loc
	if i == 1
		values[1]
	else
		percent * values[i-1] + (1 - percent) * values[i]
	end
end

# called by CMAES_RunInfo(w) during initialization
#  - produces constants stored in runInfo and accessed in runInfo.uq, runInfo.med and runInfo.lq
function locquantile(p::Float64, w::Weights)
	percent = 1 - p			# because the sort direction is "1 best" while percentile is 100% is best
	sum = 0.0
	i  = 0
	while(sum < percent)
		i += 1
		sum += w[i]
	end

	# returns the index of the percentile
	#   and the % between the that percentile and the index that holds the next better value
	#   assumes linear interpolation, consequently only produces an approximate value
	(i, (sum - percent)/w[i])
end

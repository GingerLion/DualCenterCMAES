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
    runInfo.sourceValues = SelectionSourceParms(sys) #used to be (sys)
    runInfo.allSourceValues = AllSourceParms(sys)      # for all-mirror
end

function update_scales!(runInfo::RunInfo, source::SelectionSource, state::CMAES_State)

	len = length(source.fitness)
	model = currentmodel_(state)
	μ = mu(state)
	w = Weights(normalize(map((i)->(log(μ+chrlength(state)) - log(i)), 1:μ), 1)) #less pressure than log(μ+0.5)
	#w = Weights(fill(1/μ, μ))
	orig_score = 0.0
	best_score = 0.0
	orig_scale = 0.0
	best_scale = 0.0

	for i=1:len
		(source.source[i] == :orig) ? orig_score += w[i] : best_score += w[i]
	end
	#println("sourcevalues = $(source.source)\n")
	#println("orig score = $(orig_score), best_score = $(best_score)")
	if orig_score > best_score
		orig_scale = 1 + orig_score
		best_scale = 2 - orig_scale
	elseif best_score > orig_score
		best_scale = 1 + best_score
		orig_scale = 2 - best_scale
	elseif orig_score == best_score
		orig_scale = best_scale = 1
	end

	#wrote this code because julia has some strange rounding errors e.g. 2.0000000000000000000004 == 2.0 is true but -0.0000000000000004 == 0 is fasle
	if orig_scale >= 2.0
		orig_scale = 2.0
		best_scale = 0.0
	elseif best_scale >= 2.0
		best_scale = 2.0
		orig_scale = 0.0
	end
	#println("orig_scale = $(orig_scale), best_scale = $(best_scale)")

	if (orig_scale == 2.0 || best_scale == 2.0)
		#println("54 - RunInfo_CMAES.jl: reseting the scales!")
		orig_scale!(model, 1.5)
		best_scale!(model, 0.5)
	else
		#println("orig scale = $(orig_scale), best_scale = $(best_scale)")
		orig_scale!(model, convert(Float64, orig_scale))
		best_scale!(model, convert(Float64, best_scale))
	end
	#orig_scale!(model, 1.5)
	#best_scale!(model, 0.5)
    #println("RunInfo_CMAES.jl::63 -> new orig_scale = $(orig_scale), best_scale = $(best_scale)")
	#now update SelectionSourceParms a.k.a rinfo.sourcevalues to reflect new ratio or orig_λ to best_λ
	runInfo.sourceValues = SelectionSourceParms(system(state), model)
	#println("starting gen with source = $(sourcevalues(sourcevalues(runInfo)))")
end

#=function stepsize!_(source::SelectionSource, model::CMAES_Model)
    σ_ = sigma_(model)
	#find how many of the solutions are from :best
	#do 1/5th sucess rule
    p_s = 0
	p_target = 0.5
	for i in source.source
		if i == :best p_s += 1 end
	end
	p_s = p_s / lambda(model)
	#println("p_s = $(p_s)")
	σ_ = σ_ * exp(1/3 * (p_s - p_target)/(1 - p_target))
	# set model.σ_ = σ
	sigma!_(model, σ_)
end=#
#----------------------
# RunInfo - Updates

function monitor!(nm::NoMonitor, state::CMAES_State) end

function monitor!(runInfo::RunInfo, state::CMAES_State, f::RealFitness)
  # system models
  preModel  = model(state, :pre)
  postModel = model(state, :post)
  w = weights(preModel)

  # shadowged system models
  preModel_shadow = model_(state, :pre)
  postModel_shadow = model_(state, :post)
  w_shadow = weights(preModel_shadow)

  # main system state
  parents   = population(state, :parents)
  offspring = population(state, :pre)
  selected  = population(state, :post)
  state_center = centerpopn(preModel, f)
  covMatrix = covar(postModel)

  # shadowged main system state
  parents_shadow = population_(state, :parents)
  offspring_shadow = population_(state, :pre)
  selected_shadow = population_(state, :post)
  state_center_shadow = centerpopn(preModel_shadow, f)
  state_center_shadow_ = centerpopn_(preModel_shadow, f)
  covMatrix_shadow = covar(postModel_shadow)

  # general information
  runInfo[:best]    	= best(state)
  runInfo[:covar]   	= covMatrix
  runInfo[:cond]   		= cond(covMatrix)
  runInfo[:σ]       	= sigma(state)
  runInfo[:σ_]			= sigma_(state)
  runInfo[:gen]			= gen(state)
  runInfo[:l_evals] 	= evals(state)
  try
	  (!stopEvals(state)) ? runInfo[:total_evals] = first(runInfo[:total_evals]) + evalsPerGen(system(state)) : runInfo[:total_evals] = first(runInfo[:total_evals])
  catch e
	  #if isa(e, KeyError) println("KeyError by normal") end
	  runInfo[:total_evals] = 0
  end
  #println("normal evals = $(first(runInfo[:total_evals]))")
  # general information - shadowged
  runInfo[:best_shadow] 	= best_(state)
  runInfo[:covar_shadow]   = covar_(state)
  runInfo[:cond_shadow]	= cond(covMatrix_shadow)
  runInfo[:σ_shadow]		= sigma_(state)
  runInfo[:gen_shadow]		= gen_(state)
  runInfo[:l_evals_shadow]   = evals_(state)
  try
	  (!stopEvals_(state)) ? runInfo[:total_evals_] = first(runInfo[:total_evals_]) + evalsPerGen_(system(state)) : runInfo[:total_evals_] = first(runInfo[:total_evals_])
  catch e
	  #if isa(e, KeyError)	println("KeyError by dualcenter") end
	  runInfo[:total_evals_] = 0
  end
  #println("dualcenter evals = $(first(runInfo[:total_evals_]))")
  # system state information
  runInfo[:center]  	= (state_center[:chr, 1], state_center[:fit, 1])
  runInfo[:fitsummary]	= fitsummary(selected, runInfo, w)
  src = SelectionSource(sourcevalues(runInfo), state)
  runInfo[:source] = deepcopy(src)
  runInfo[:center_shadow_source] = first(src.source)

  update_scales!(runInfo, src, state)

  #stepsize!_(src, postModel_shadow)
  #online update of model parms
  #update!(src, postModel, system(state))
  # system state information - shadowged
  runInfo[:best_fit] = bestfitness(state)
  runInfo[:best_shadow_fit] = bestfitness_(state)
  runInfo[:center_shadow]  = (state_center_shadow[:chr, 1], state_center_shadow[:fit, 1])

  runInfo[:firstRequest] = firstRequest(state)
  #   Mirroring state
  #   note: currently does not include mirror of covar, signma or paths,
  #		    i.e. only center is mirrored from the model
  #   note: different calculation for center than state-center

  # offspring only mirror
  (offspring_selected, offs_sortOrder) = truncation(offspring, mu(state))
  offspring_center = center(offspring_selected, w)
  runInfo[:center_fit] = state_center[:fit,1]
  runInfo[:center_fit_shadow] = state_center_shadow[:fit,1]
  runInfo[:center_fit_shadow_] = state_center_shadow_[:fit,1]
  #println("center_ = $(state_center_shadow_[:fit,1])")
  #println("center = $(runInfo[:center_fit_shadow]), center_ = $(runInfo[:center_fit_shadow_])")
  runInfo[:offspring_center] = (offspring_center, objfn(f)(offspring_center))
  runInfo[:offs_fitsummary]	 = fitsummary(offspring_selected,runInfo,  w)

  # all (parents + center + offspring) mirror
  #(all_selected, all_sortOrder) = truncation(pcat(parents, state_center, offspring), mu(state))
  #all_center = center(all_selected, w)

  #runInfo[:all_center]		 = (all_center, objfn(f)(all_center))
  #runInfo[:all_fitsummary]   = fitsummary(all_selected, runInfo, w)
  #runInfo[:all_source]       = SelectionSource(allsourcevalues(runInfo), all_sortOrder)

  #eigenvalue & eigenvector study
  if !eigenerror(state) && !eigenerror_(state)
	  #NON-shadowGED
	  eigvals = eigendecomp(postModel)[1]
	  eigvecs = eigendecomp(postModel)[2]
	  #runInfo[:Dinv] = diag(inv(diagm(0 => eigvals)),0)
	  runInfo[:eig_vals] = eigvals
	  runInfo[:eig_vecs] = eigvecs
	  runInfo[:BD] = *(eigvecs,eigvals)
	  runInfo[:eig_min] = minimum(eigvals)
	  runInfo[:eig_max] = maximum(eigvals)
	  #shadowGED
	  eigvals_shadow = eigendecomp(postModel_shadow)[1]
	  eigvecs_shadow = eigendecomp(postModel_shadow)[2]
	  #runInfo[:Dinv_shadow] = diag(inv(diagm(0 => eigvals_shadow)),0)
	  runInfo[:eig_vals_shadow] = eigvals_shadow
	  runInfo[:eig_vecs_shadow] = eigvecs_shadow
	  runInfo[:BD_shadow] = *(eigvecs_shadow,eigvals_shadow)
	  runInfo[:eig_min_shadow] = minimum(eigvals_shadow)
	  runInfo[:eig_max_shadow] = maximum(eigvals_shadow)
  else
	  #NON-shadowGED
	  runInfo[:eig_vals] = fill(NaN,chrlength(postModel))
	  runInfo[:eig_vecs] = fill(NaN,(chrlength(postModel),chrlength(postModel)))
	  #runInfo[:Dinv] = fill(NaN,chrlength(postModel))
	  runInfo[:BD] = fill(NaN,(chrlength(postModel)))
	  runInfo[:eig_min] = NaN
	  runInfo[:eig_max] = NaN
	  #shadowGED
	  runInfo[:eig_vals_shadow] = fill(NaN,chrlength(postModel_shadow))
	  runInfo[:eig_vecs_shadow] = fill(NaN,(chrlength(postModel_shadow),chrlength(postModel_shadow)))
	  #runInfo[:Dinv] = fill(NaN,chrlength(postModel))
	  runInfo[:BD_shadow] = fill(NaN,(chrlength(postModel_shadow)))
	  runInfo[:eig_min_shadow] = NaN
	  runInfo[:eig_max_shadow] = NaN

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
	#shadowged model parms
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

	#paths - shadowged
	runInfo[:p_σ_shadow] = p_σ(postModel_shadow)
	runInfo[:p_c_shadow] = p_c(postModel_shadow)
	runInfo[:y_shadow] = weightedavg(state.sW_shadow)
	runInfo[:invsqrtC_shadow] = invsqrtC(postModel_shadow)
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

	#update equations in parts - shadowged
	runInfo[:p_σ_part1_shadow] = p_σ_part1(postModel_shadow)
	runInfo[:p_σ_part2_shadow] = p_σ_part2(postModel_shadow)
	#runInfo[:p_σ_part3_shadow] = p_σ_part3(postModel_shadow,state.sW_shadow)

	runInfo[:σ_part1_shadow] = σ_part1(postModel_shadow)
	runInfo[:σ_part2_shadow] = σ_part2(postModel_shadow)
	runInfo[:σ_part3_shadow] = σ_part3(postModel_shadow)
	#runInfo[:σ_part4_shadow] = σ_part4(postModel_shadow)
	#runInfo[:σ_part5_shadow] = σ_part5(postModel_shadow)
	runInfo[:σ_part6_shadow] = σ_part6(postModel_shadow)

	runInfo[:p_c_part1_shadow] = p_c_part1(postModel_shadow)
	runInfo[:p_c_part2_shadow] = p_c_part2(postModel_shadow)
	runInfo[:p_c_part3_shadow] = p_c_part3(postModel_shadow,state.sW_shadow)

	runInfo[:C_part1_shadow] = C_part1(postModel_shadow)
	runInfo[:C_part2_shadow] = C_part2(postModel_shadow)
	runInfo[:C_part3_shadow] = C_part3(postModel_shadow,state.sW_shadow)
end

function monitor!(runInfo::RunInfo, state::CMAES_State, restart::RestartState, f::RealFitness)
  monitor!(runInfo, state, f)
  runInfo[:g_evals] = evals(restart)
  runInfo[:restart]	= rep(restart)
  runInfo[:lambda]  = lambda(state)
  #runInfo[:lambda_shadow] = lambda_(state)
end

#-----------------------
# Helper functions

center(popn::SortedPopulation, w::Weights) = vec(mean(popn[:chr, :], w, dims=2))
centerpopn(model::CMAES_Model, f::RealFitness) = RegularPopulation(centermember(model, objfn(f)))
centerpopn_(model::CMAES_Model, f::RealFitness) = RegularPopulation(centermember_(model, objfn(f)))

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

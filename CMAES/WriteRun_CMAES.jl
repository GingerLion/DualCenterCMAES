#----------------------
# Writing RunInfo

function write_run(returnInfo::ReReturnInfo, sys::CMAES_System, fitFn::RealFitness;
                   prefixValues = [], prefixNames = [], initialize = false,
                   fileName = "ipoprun", path = "", sep = "\t", ext = ".txt")
  ext = (sep == ",") ? ".csv" : ext
  fileName = "$(path)/$fileName$ext"
  if initialize
    write_run_header(fileName, returnInfo, prefixNames, sep)
  end
  write_run(fileName, returnInfo, sys, fitFn, prefixValues, sep)
end

function write_run_header(fileName::String, returnInfo::ReReturnInfo, prefixNames::Vector{String}, sep::String)
  open(fileName, "w") do f
    for name in prefixNames
      write(f, "$name$sep")
	end
    write(f, "restart$(sep)")
    write(f, "lambda$(sep)")
    #write(f, "lambda_shadow$(sep)")
    write(f, "gen$(sep)")
    #write(f, "gen_shadow$(sep)")
    write(f, "evals$(sep)")
    write(f, "cond$(sep)")
    write(f, "cond_shadow$(sep)")
    write(f, "sigma$(sep)")
    write(f, "sigma_shadow$(sep)")
    #write(f, "eig_min$(sep)")
    #write(f, "eig_max$(sep)")
    #write(f, "sys_worstGene$(sep)")
    #write(f, "sys_euclDist$(sep)")
    #write(f, "sys_bestGene$(sep)")
    #write(f, "sys_centerFit$(sep)")
    #write(f, "sys_maxFit$(sep)")
    #write(f, "sys_wUQFit$(sep)")
    #write(f, "sys_wMedFit$(sep)")
    #write(f, "sys_wAvgFit$(sep)")
    #write(f, "sys_wLQFit$(sep)")
    #write(f, "sys_MinFit$(sep)")
    write(f, "sys_p$(sep)")
    write(f, "sys_c$(sep)")
    write(f, "sys_o$(sep)")
    #write(f, "sys_pMaxRank$(sep)")
    #write(f, "sys_pAvgRank$(sep)")
    #write(f, "sys_pMinRank$(sep)")
    #write(f, "sys_ctrRank$(sep)")
    #write(f, "sys_oMaxRank$(sep)")
    #write(f, "sys_oAvgRank$(sep)")
    #write(f, "sys_oMinRank$(sep)")
    #write(f, "sys_pMaxFit$(sep)")
    #write(f, "sys_pAvgFit$(sep)")
    #write(f, "sys_pMinFit$(sep)")
    #write(f, "sys_ctrFit$(sep)")
    #write(f, "sys_oMaxFit$(sep)")
    #write(f, "sys_oAvgFit$(sep)")
    #write(f, "sys_oMinFit$(sep)")
    #write(f, "offsp_worstGene$(sep)")
    #write(f, "offsp_euclDist$(sep)")
    #write(f, "offsp_bestGene$(sep)")
    #write(f, "offsp_centerFit$(sep)")
    #write(f, "offsp_maxFit$(sep)")
    #write(f, "offsp_wUQFit$(sep)")
    #write(f, "offsp_wMedFit$(sep)")
    #write(f, "offsp_wAvgFit$(sep)")
    #write(f, "offsp_wLQFit$(sep)")
    #write(f, "offsp_MinFit$(sep)")
    #write(f, "all_worstGene$(sep)")
    #write(f, "all_euclDist$(sep)")
    #write(f, "all_bestGene$(sep)")
    #write(f, "all_centerFit$(sep)")
    #write(f, "all_maxFit$(sep)")
    #write(f, "all_wUQFit$(sep)")
    #write(f, "all_wMedFit$(sep)")
    #write(f, "all_wAvgFit$(sep)")
    #write(f, "all_wLQFit$(sep)")
    #write(f, "all_MinFit$(sep)")
    #write(f, "all_p$(sep)")
    #write(f, "all_c$(sep)")
    #write(f, "all_o$(sep)")
    #write(f, "all_pMaxRank$(sep)")
    #write(f, "all_pAvgRank$(sep)")
    #write(f, "all_pMinRank$(sep)")
    #write(f, "all_ctrRank$(sep)")
    #write(f, "all_oMaxRank$(sep)")
    #write(f, "all_oAvgRank$(sep)")
    #write(f, "all_oMinRank$(sep)")
    #write(f, "all_pMaxFit$(sep)")
    #write(f, "all_pAvgFit$(sep)")
    #write(f, "all_pMinFit$(sep)")
    #write(f, "all_ctrFit$(sep)")
    #write(f, "all_oMaxFit$(sep)")
    #write(f, "all_oAvgFit$(sep)")
    #write(f, "all_oMinFit$(sep)")
    #write(f, "center$(sep)")
    #write(f, "eig_vecs\n")
    #non-vectors
    #parameters
    write(f, "c_sig$(sep)")
    write(f, "d_sig$(sep)")
    write(f, "mu_eff$(sep)")
    write(f, "c_1$(sep)")
    write(f, "c_mu$(sep)")
    write(f, "c_c$(sep)")
    write(f, "h_sig$(sep)")
    #stepsize update equations
    write(f, "c_sig/d_sig$(sep)")
    write(f, "||p_sig||$(sep)")
    write(f, "c_sig/d_sig_shadow$(sep)")
    write(f, "||p_sig||_shadow$(sep)")
    #write(f, "chi_mean$(sep)")
    #write(f, "exp(c_sig/d_sig(||p_sig|| / chi_mean) -1)$(sep)")
    #write(f, "eig_vals$(sep)\n")
    #vectors

    write(f, "best_fit$(sep)")
    write(f, "best_fit_shadow$(sep)")
    write(f, "centerFit$(sep)")
    write(f, "centerFitshadow$(sep)")
    write(f, "center$(sep)")
    write(f, "center_shadow$(sep)")
    #write(f, "B$(sep)")
    #write(f, "D$(sep)")
    write(f, "BD$(sep)")
    write(f, "BD_shadow$(sep)")
    #sigma path update equations
    #write(f, "BD^-1B^T$(sep)")
    write(f, "y$(sep)")
    write(f, "y_shadow$(sep)")
    #write(f, "(1-c_sig)p_sig$(sep)")
    #write(f, "sqrt(c_sig(2 - c_sig)mu_eff)$(sep)")
    #write(f, "sqrt(c_sig(2 - c_sig)mu_eff)BD^-1B^Ty$(sep)")
    write(f, "p_sig$(sep)")
    write(f, "p_sig_shadow$(sep)")
    #cov path update equations
    #write(f, "(1-c_c)p_c$(sep)")
    #write(f, "h_sig(sqrt(c_c(2-c_c)mu_eff))$(sep)")
    #write(f, "h_sig(sqrt(c_c(2-c_c)mu_eff))y$(sep)")
    write(f, "p_c$(sep)")
    write(f, "p_c_shadow$(sep)")
    #cov matrix update equations
    #write(f, "(1-c_1-c_mu)C$(sep)")
    #write(f, "c_1(p_c)p_c'$(sep)")
    #write(f, "c_mu + mu-update$(sep)")
    write(f, "C$(sep)")
    write(f, "C_shadow\n")
  end
end

function write_run(fileName::String, returnInfo::ReReturnInfo, sys::CMAES_System,
                   fitFn::RealFitness, prefixValues::Vector, sep::String)
  ri = returnInfo.runInfo
  open(fileName, "a") do f
    for i = 1:(length(ri)-1)
  	  for value in prefixValues
  	    write(f, "$value$sep")
  	  end
      write(f, "$(ri[:restart][i])$sep")
      write(f, "$(ri[:lambda][i])$sep")
      #write(f, "$(ri[:lambda_shadow][i])$sep")
      write(f, "$(ri[:gen][i])$sep")
      #write(f, "$(ri[:gen_][i])$sep")
      write(f, "$(ri[:g_evals][i])$sep")
      write(f, "$(ri[:cond][i])$sep")
      write(f, "$(ri[:cond_shadow][i])$sep")
      write(f, "$(ri[:σ][i])$sep")
      write(f, "$(ri[:σ_shadow][i])$sep")
      #write(f, "$(ri[:eig_min][i])$sep")
      #write(f, "$(ri[:eig_max][i])$sep")
      #write(f, "$(ri[:eig_min_shadow][i])$sep")
      #write(f, "$(ri[:eig_max_shadow][i])$sep")
      #write_centersummary(f, ri[:center][i], fitFn, sep)
      #write_fitsummary(f, ri[:fitsummary][i], sep)
      write_sourcesummary(f, ri[:source][i], sep)
      #write_centersummary(f, ri[:offspring_center][i], fitFn, sep)
      #write_fitsummary(f, ri[:offs_fitsummary][i], sep)
      #write_centersummary(f, ri[:all_center][i], fitFn, sep)
      #write_fitsummary(f, ri[:all_fitsummary][i], sep)
      #write_sourcesummary(f, ri[:all_source][i], sep)
      #write(f, "$(ri[:eig_vals][i])$sep")
      #write(f, "$(ri[:eig_vecs][i])\n")
      #write(f, "$(ri[:covar][i])\n")
      write(f, "$(ri[:c_σ][i])$sep")
      write(f, "$(ri[:d_σ][i])$sep")
      write(f, "$(ri[:μ_eff][i])$sep")
      write(f, "$(ri[:c_1][i])$sep")
      write(f, "$(ri[:c_μ][i])$sep")
      write(f, "$(ri[:c_c][i])$sep")
      write(f, "$(ri[:h_σ][i])$sep")

      write(f, "$(ri[:σ_part1][i])$sep")
      write(f, "$(ri[:σ_part2][i])$sep")
      write(f, "$(ri[:σ_part1_shadow][i])$sep")
      write(f, "$(ri[:σ_part2_shadow][i])$sep")
      #write(f, "$(ri[:σ_part3][i])$sep")
      #write(f, "$(ri[:σ_part4][i])$sep")
      #write(f, "$(ri[:σ_part5][i])$sep")
     # write(f, "$(ri[:σ_part6][i])$sep")
      #write(f, "$(ri[:Dinv][i])$sep")
      write(f, "$(ri[:best_fit][i])$sep")
      write(f, "$(ri[:best_shadow_fit][i])$sep")
      write(f, "$(ri[:center_fit][i])$sep")
      write(f, "$(ri[:center_fit_shadow][i])$sep")
      write(f, "$(ri[:center][i][1])$sep")
      write(f, "$(ri[:center_shadow][i][1])$sep")
      #write(f, "$(ri[:eig_vecs][i])$sep")
      #write(f, "$(ri[:eig_vals][i])$sep")
      write(f, "$(ri[:BD][i])$sep")
      write(f, "$(ri[:BD_shadow][i])$sep")
      #write(f, "$(ri[:invsqrtC][i])$sep")
      write(f, "$(ri[:y][i])$sep")
      write(f, "$(ri[:y_shadow][i])$sep")

      #write(f, "$(ri[:p_σ_part1][i])$sep")
      #write(f, "$(ri[:p_σ_part2][i])$sep")
      #write(f, "$(ri[:p_σ_part3][i])$sep")
      write(f, "$(ri[:p_σ][i])$sep")
      write(f, "$(ri[:p_σ_shadow][i])$sep")

      #write(f, "$(ri[:p_c_part1][i])$sep")
      #write(f, "$(ri[:p_c_part2][i])$sep")
      #write(f, "$(ri[:p_c_part3][i])$sep")
      write(f, "$(ri[:p_c][i])$sep")
      write(f, "$(ri[:p_c_shadow][i])$sep")
      #write(f, "$(ri[:C_part1][i])$sep")
      #write(f, "$(ri[:C_part2][i])$sep")
      #write(f, "$(ri[:C_part3][i])$sep")
      write(f, "$(ri[:covar][i])$sep")
      write(f, "$(ri[:covar_shadow][i])\n")
      #write(f, "$(ri[:eig_vals][i])\n")
    end
  end
end

function write_sourcesummary(f::IO, source::SelectionSource, sep::String)
	write_proportions(f, source, sep)
	#write_rankmarginals(f, source, sep)
	#write_fitmarginals(f, source, sep)
end

function write_proportions(f::IO, source::SelectionSource, sep::String)
  (parent, center, offspring) = proportions(source)
  write(f, "$parent$sep")
  write(f, "$center$sep")
  write(f, "$offspring$sep")
end

write_rankmarginals(f::IO, source::SelectionSource, sep::String) = write_sourcesummary(f, source, rankmarginals, sep)
write_fitmarginals(f::IO, source::SelectionSource, sep::String) = write_sourcesummary(f, source, fitmarginals, sep)

function write_sourcesummary(f::IO, source::SelectionSource, summaryfn::Function, sep::String)
  (pmax, pavg, pmin, center, omax, oavg, omin) = summaryfn(source)
  write(f, "$pmax$sep")
  write(f, "$pavg$sep")
  write(f, "$pmin$sep")
  write(f, "$center$sep")
  write(f, "$omax$sep")
  write(f, "$oavg$sep")
  write(f, "$omin$sep")
end

function write_fitsummary(f::IO, fitValues::Tuple, sep::String)
  (max, uq, med, avg, lq, min) = fitValues
  write(f, "$max$sep")
  write(f, "$uq$sep")
  write(f, "$med$sep")
  write(f, "$avg$sep")
  write(f, "$lq$sep")
  write(f, "$min$sep")
end

function write_centersummary(f::IO, center::Tuple, fitFn::RealFitness, sep::String)
  centerDiff = center[1] - optimum(fitFn)
  fit = center[2]
  euclDist = norm(centerDiff)
  worstGene = norm(centerDiff, Inf)
  bestGene = norm(centerDiff, -Inf)
  write(f, "$worstGene$sep")
  write(f, "$euclDist$sep")
  write(f, "$bestGene$sep")
  write(f, "$fit$sep")
end

#---------------------

function write_selectionsource(sys::CMAES_System, ri::RunInfo; fileName = "cmaesselectionsource",
	                           path = "", sep = "\t", ext = ".txt")
  ext = (sep == ",") ? ".csv" : ext
  open("$(path)/$fileName$ext", "w") do f
    write(f, "gen$(sep)evals$(sep)sigma$(sep)")
    write(f, "p%$(sep)c%$(sep)o%$(sep)")
    write(f, "p_rank$(sep)c_rank$(sep)o_rank$(sep)")
    write(f, "p_fit$(sep)c_fit$(sep)o_fit\n")
    for i = 1:(length(ri)-1)
      write(f, "$(ri[:gen][i])$(sep)$(ri[:l_evals][i])$(sep)$(ri[:σ][i])$(sep)")
      (pPercent, cPercent, oPercent) = proportions(sys, ri[:source][i])
      write(f, "$pPercent$(sep)$cPercent$(sep)$oPercent$(sep)")
      (pAvgRank, cRank, oAvgRank) = avgrank(ri[:source][i])
      write(f, "$pAvgRank$(sep)$cRank$(sep)$oAvgRank$(sep)")
      (pAvgFit, cFit, oAvgFit) = avgfit(ri[:source][i])
      write(f, "$pAvgFit$(sep)$cFit$(sep)$oAvgFit$(sep)")
      write(f, "$(ri[:center][i])")
      write(f, "\n")
    end
  end
end

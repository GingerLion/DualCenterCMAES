#----------------------------------------------------------
#  runInfo::Monitor
#----------------------------------------------------------

function update!(nm::NoMonitor, sys::DE_System) end

function monitor!(runInfo::Monitor, state::DE_State, fit::Fitness)
  runInfo[:popn]  = deepcopy(state.popn)
  runInfo[:targets] = deepcopy(state.contenders.targets)
  runInfo[:samples]  = deepcopy(state.contenders.samples)
  runInfo[:alternates] = deepcopy(state.contenders.alternates)
  map(runInfo[:alternates]) do alternate 
    evaluate!(alternate, fit.objFn, fit.direction)
  end
end

function monitor!(runInfo::RunInfo, state::DE_State, restart::RestartState, f::Fitness)
  monitor!(runInfo, state, f)
  runInfo[:g_evals] = evals(restart)
  runInfo[:rep]   = rep(restart)
  runInfo[:popn_size]  = lambda(state)
end

# function write_selectionsource(sys::CMAES_System, rsys::CMAES_Restart, ri::RunInfo; fileName = "ipopselectionsource", path = "", sep = "\t", ext = ".txt")
#   ext = (sep == ",") ? ".csv" : ext
#   open("$(path)/$fileName$ext", "w") do f
#     write(f, "rep$(sep)popnSize$(sep)gen$(sep)evals$(sep)")
#     write(f, "p%$(sep)c%$(sep)o%$(sep)")
#     write(f, "p_rank$(sep)c_rank$(sep)o_rank$(sep)")
#     write(f, "p_fit$(sep)c_fit$(sep)o_fit\n")
#     for i = 1:(length(ri)-1)
#       write(f, "$(ri[:rep][i])$(sep)$(ri[:popn_size][i])$(sep)$(ri[:gen][i])$(sep)$(ri[:g_evals][i])$(sep)")
#       (pPercent, cPercent, oPercent) = proportions(sys, ri[:source][i])
#       write(f, "$pPercent$(sep)$cPercent$(sep)$oPercent$(sep)")
#       (pAvgRank, cRank, oAvgRank) = avgrank(ri[:source][i])
#       write(f, "$pAvgRank$(sep)$cRank$(sep)$oAvgRank$(sep)")
#       (pAvgFit, cFit, oAvgFit) = avgfit(ri[:source][i])
#       write(f, "$pAvgFit$(sep)$cFit$(sep)$oAvgFit$(sep)")
#       write(f, "$(ri[:center][i])")
#       write(f, "\n")
#     end
#   end
# end

function write_selectionsource(sys::DE_System, ri::RunInfo; fileName = "deselectionsource", path = "", sep = "\t", ext = ".txt")
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
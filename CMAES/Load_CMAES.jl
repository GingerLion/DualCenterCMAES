base_path = "/Users/dillon/v0.4.06"
cmaes_path = "$(base_path)/CMAES"

if !@isdefined commonSourceLoaded
  include("$(base_path)/Common/Load_Common.jl")
end

cmaes_files = ["Types_CMAES.jl",
               "Verbose_CMAES.jl",
               "Model_CMAES.jl",
               "Noise_CMAES.jl",
               "Selection_CMAES.jl",
               "SelectionSource_CMAES.jl",
               "System_CMAES.jl",
               "State_CMAES.jl",
               "Evolve_CMAES.jl",
               "Restart_CMAES.jl",
               "RestartCriteria_CMAES.jl",
               "RunInfo_CMAES.jl",
               "WriteRun_CMAES.jl"]

load_files(cmaes_files; path = cmaes_path)

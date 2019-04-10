base_path = "/Users/mwineberg/Dropbox/2 Work/Research/Projects/03 Programs/EAs in Julia/Production/v0.4.02"
de_path = "$(base_path)/DE"

if !@isdefined commonSourceLoaded
  include("$(base_path)/Common/Load_Common.jl")
end

de_files = ["Types_DE.jl",
            "Verbose_DE.jl",
            "Initializers_DE.jl",
            "Samplers_DE.jl",
            "Reproduction_DE.jl",
            "Selection_DE.jl",
            "System_DE.jl",
            "State_DE.jl",
            "Evolve_DE.jl",
            "Restart_DE.jl",
            "RunInfo_DE.jl"]


load_files(de_files; path = de_path)

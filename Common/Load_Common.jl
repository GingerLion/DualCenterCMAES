using Pkg
pkg"up"
pkg"add Distributions"
pkg"add DataStructures"
pkg"add Distances"
pkg"add Polynomials"
pkg"add StatsBase"

using DataStructures, Distributions, Distances, StatsBase, Polynomials, Statistics, LinearAlgebra
import Base.getindex, Base.setindex!, Base.+, Base.*
import Base.print, Base.println, Base.length, Base.size
import Base.argmin, Base.argmax, Base.maximum, Base.minimum, Base.count
import Statistics.median, Statistics.quantile
import Distances.evaluate, DataStructures.update!
import StatsBase.weights

function load_files(files; path="")
	incurrent() = (path == "")
	includefile(file) = include(incurrent() ? file : "$(path)/$file")

	for file in files
		includefile(file)
	end
end
base_path = "C:/Users/dillo/OneDrive - University of Guelph/Research/EA_Julia_Wineberg/v0.4.06"
common_path = "$(base_path)/Common"

common_files = ["EA_Types.jl",
				"EA_Verbose.jl",
            	"EA_Fitness.jl",
            	"EA_BestFitHistory.jl",
            	"EA_Population.jl",
            	"EA_Population#Regular.jl",
            	"EA_Population#Sorted.jl",
            	"EA_Population#Member.jl",
            	"EA_Selection.jl",
            	"EA_Placement.jl",
            	"EA_RestartState.jl",
            	"EA_State.jl",
            	"EA_Evolve.jl",
            	"EA_TestFunctions.jl",
            	"EA_ReturnInfo.jl",
            	"EA_ReturnInfoWrite.jl",
            	"EA_RunEA.jl",
            	"EA_Monitoring.jl"]

load_files(common_files; path = common_path)
commonSourceLoaded = true

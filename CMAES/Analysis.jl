using Pkg
Pkg.update()
using CSV, DataFrames, Statistics
using Plots

# file path: "C:/Users/dillo/OneDrive - University of Guelph/Research/EA_Julia_Wineberg/v0.4.02/Experiments/"

#fn: eigval_analysis - plots a graph of the mean eigenvalue vector
#params :
# fileName::String - full file path; MUST be '/t' delimited file (for df manipulation reasons)
# dim::Int64 - dimension of chromosomes
# elitism::Bool - apply analysis for elitism being true or false (because this is the goal of our research)
# duringStall::Bool keyword arg - if set to true, does analysis when CMAES is stalling e.g. (cond > 2000, sigma > 1000, all_p or sys_p == 1)
function mean_eigenvalues(fileName::String,fnName::String,dim::Int64, elitism::Bool; duringStall = false, cond = 2000, sigma = 500)
    #read the data in from the tab delimited file into 'df' as a DataFrame
    df = CSV.read(fileName,delim = '\t',header=true, use_mmap=true) #use_mmap=true for speedup on Windows OS
    #based on parameters extract data from the DataFrame so the operations below aren't performed on irrelevant data
    if(duringStall & elitism)
        targDF = df[(df.fn .== fnName) .& (df.dim .== dim) .& (df.elitism .== elitism) .& (df.cond .> cond) .& (df.sigma .> sigma) .& (df.sys_p .== 1) .& (map(x -> !(isnan(x)), df.eig_min)),:]
    elseif (duringStall)
        targDF = df[(df.fn .== fnName) .& (df.dim .== dim) .& (df.elitism .== elitism) .& (df.cond .> cond) .& (df.sigma .> sigma) .& (map(x -> !(isnan(x)), df.eig_min)),:]
    else
        targDF = df[(df.fn .== fnName) .& (df.dim .== dim) .& (df.elitism .== elitism) .& (map(x -> !(isnan(x)),df.eig_min)),:]
    end

    #remove "[]" brackets from the eigenvalue vectors which for some reason reads into the DataFrame as a String
    targDF[:eig_vals] = map(x -> replace(x,"[" => ""), targDF[:eig_vals])
    targDF[:eig_vals] = map(x -> replace(x,"]" => ""), targDF[:eig_vals])
    #remove ","'s from the eig_val arrays and convert them to Array{Array{Float64,1},1}
    targDF[:eig_vals] = map(x -> [parse(Float64,ss) for ss in split(x,",")],targDF[:eig_vals])

    #you cannot perform mean operations unless the dataset is Array{Float,2} so......
    #extract eig_vals
    eigvals = hcat(targDF[:eig_vals]...)'
    mean_eigvals = mean(eigvals, dims=1)
    #mean_eigvals reduces to a single vector but still is still of type Array{Float,2}
    #convert mean_eigvals to Array{Float64,1} for plotting
    temp = mean_eigvals
    mean_eigvals = Array{Float64,1}(undef,dim)
    for i=1:size(temp,2)
        mean_eigvals[i] = temp[i]
    end
    #plot the mean eigenvalues
    x = collect(1.0:convert(Float64,dim)) #x-axis with a point for each dimension
    plotly() #use this plotting backend
    #plot graph
    plot(x,mean_eigvals,title="Mean Eigenvalues of [$fnName $dim elitism=$elitism duringStall=$duringStall cond > $cond  sigma > $sigma]",label="mean eigenvalues",xlabel="index",ylabel="eigenvalue")
end

# fileName = "C:/Users/dillo/OneDrive - University of Guelph/Research/EA_Julia_Wineberg/v0.4.05/Experiments/ipop+_run#allfns_crn_shadow.txt"

#fn: rav analysis
#params :
# fileName::String - full file path; MUST be '/t' delimited file (for df manipulation reasons)
# dim::Int64 - dimension of chromosomes
# elitism::Bool - apply analysis for elitism being true or false (because this is the goal of our research)
# restart::Int64 - the restart run to apply the analysis to (0-?)
function rav_analysis(fileName::String, fn::String, dim::Int64)
    #read the data in from the tab delimited file into 'df' as a DataFrame
    df = CSV.read(fileName,delim = '\t',header=true, use_mmap=true) #use_mmap=true for speedup on Windows OS

    targDF = df[(df.fn .== fn) .& (df.dim .== dim) .& (df.elitism .== true),:]

    rav_collection = fill(Float64,size(targDF))

    #remove "[]" brackets from the center & noise vectors which for some reason reads into the DataFrame as a String
    for symbol in [:center, :center_shadow, :y, :y_shadow]
        targDF[symbol] = map(x -> replace(x,"[" => ""), targDF[symbol])
        targDF[symbol] = map(x -> replace(x,"]" => ""), targDF[symbol])
        #remove ","'s from the arrays and convert them to Array{Array{Float64,1},1}
        targDF[symbol] = map(x -> [parse(Float64,ss) for ss in split(x,",")],targDF[symbol])
    end
    #use dictionary to store the Arrays -since i have no idea how to convert Array{Array{Float64,1},1} to Array{Float64,2}
    ctr_storage = Dict{Int64,Array{Float64,1}}()
    ctr_shadow_storage = Dict{Int64,Array{Float64,1}}()
    y_storage = Dict{Int64,Array{Float64,1}}()
    y_shadow_storage = Dict{Int64,Array{Float64,1}}()
    for i = 1:size(targDF[:center],1)
        ctr_storage[i] = targDF[:center][i]
        ctr_shadow_storage[i] = targDF[:center_shadow][i]
        y_storage[i] = targDF[:y][i]
        y_shadow_storage[i] = targDF[:y_shadow][i]
    end
    #calculate ravs
    #declare Dict to hold ravs
    rav_storage = Dict{Int64,Array{Float64,1}}()
    rav_storage[1] = fill(NaN,size(targDF[:center],1))
    rav_shadow_storage = Dict{Int64,Array{Float64,1}}()
    rav_shadow_storage[1] = fill(NaN,size(targDF[:center],1))
    #arrays to hold norms of ravs
    norm_rav = fill(0.0,length(keys(ctr_storage)))
    norm_rav_shadow = fill(0.0,length(keys(ctr_shadow_storage)))
    #arrays to count the # of negative and positive values in each rav
    num_of_plus = zeros(length(keys(ctr_storage)))
    num_of_minus = zeros(length(keys(ctr_storage)))
    num_of_plus_shadow = zeros(length(keys(ctr_storage)))
    num_of_minus_shadow = zeros(length(keys(ctr_storage)))

    for i in 2:length(keys(ctr_storage))
        rav_storage[i] = (ctr_storage[i]) - (ctr_storage[i-1] + (targDF[:sigma][i-1] * y_storage[i-1]))
        rav_shadow_storage[i] = (ctr_shadow_storage[i]) - (ctr_shadow_storage[i-1] + (targDF[:sigma_shadow][i-1] * y_shadow_storage[i-1]))
        norm_rav[i] = norm(rav_storage[i])
        norm_rav_shadow[i] = norm(rav_shadow_storage[i])
    end

    for i in 1:length(keys(ctr_storage))
        for j in 1:dim
            rav_storage[i][j] < 0  ?   num_of_minus[i] += 1 :  num_of_plus[i] += 1
            rav_shadow_storage[i][j] < 0  ?   num_of_minus_shadow[i] += 1 :  num_of_plus_shadow[i] += 1
        end
    end

    #write this output to csv file
    open("$fn-rav_analysis.csv", "w") do f
        write(f, "fn,")
        write(f, "norm_rav,")
        write(f, "norm_rav_shadow,")
        write(f, "num_plus,")
        write(f, "num_minus,")
        write(f, "num_plus_shadow,")
        write(f, "num_minus_shadow,")
        write(f, "RAV,")
        write(f, "RAV_shadow\n")
        for i = 1:length(keys(rav_storage))
            write(f, "$(fn),")
            write(f, "$(norm_rav[i]),")
            write(f, "$(norm_rav_shadow[i]),")
            write(f, "$(num_of_plus[i]),")
            write(f, "$(num_of_minus[i]),")
            write(f, "$(num_of_plus_shadow[i]),")
            write(f, "$(num_of_minus_shadow[i]),")
            write(f, "$(rav_storage[i]),")
            write(f, "$(rav_shadow_storage[i])\n")
        end
    end
end
fileName = "C:/Users/dillo/OneDrive - University of Guelph/Research/EA_Julia_Wineberg/v0.4.05/Experiments/ipop+_run#allfns_crn_shadow.txt"
function rav_expr(fileName::String, dim)
    for fn in ["rastrigin","levy","elliptical","ackley","griewank"]
        rav_analysis(fileName, fn, dim)
    end
end

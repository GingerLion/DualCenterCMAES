# struct SelectonSourceParms
#     sourceValues::Vector{Symbol}
#     selectionType::Type
# end

selectiontype(sParms::SelectionSourceParms) = sParms.selectionType
allsources(sParms::SelectionSourceParms) = sParms.sourcevalues          # depreciated
sourcevalues(sParms::SelectionSourceParms) = sParms.sourceValues

function SelectionSourceParms(sys::CMAES_System)
    SelectionSourceParms(elitism(sys), includecenter(sys), mu(sys), lambda(sys), η(sys))
end

function SelectionSourceParms(sys::CMAES_System, model::CMAES_Model)
    SelectionSourceParms(orig_scale(model), best_scale(model), elitism(sys), includecenter(sys), mu(sys), lambda(sys), η(sys))
end

function AllSourceParms(sys::CMAES_System)
    λ_1 = floor(Int64, lambda(sys)/2)
    λ_2 = lambda(sys) - λ_1
    AllSourceParms(mu(sys), λ_1, λ_2, η(sys))
end
EliteSourceParms(sys::CMAES_System) = EliteSourceParms(mu(sys), lambda(sys), η(sys))
CenterSelSourceParms(sys::CMAES_System) = CenterSelSourceParms(lambda(sys))

function SelectionSourceParms(elitism::Bool, includeCenter::Bool, μ::Int, λ::Int, η::Int)
  first_half = floor(Int64, λ/2)
  second_half = λ - first_half
  if elitism && includeCenter
    AllSourceParms(μ, first_half, second_half, η)
  elseif elitism
    EliteSourceParms(μ, first_half, second_half, η)
  elseif includeCenter
    CenterSelSourceParms(first_half, second_half)
  else
    SingleSourceParms(first_half, second_half)
  end
end

function SelectionSourceParms(orig_scale::Float64, best_scale::Float64, elitism::Bool, includeCenter::Bool, μ::Int, λ::Int, η::Int)
    if orig_scale == 2.0
        first_half = λ
        second_half = 0
    else
        first_half = floor(Int64, orig_scale * floor(Int64,λ/2))
        second_half = λ - first_half
    end
  if elitism && includeCenter
    AllSourceParms(μ, first_half, second_half, η)
  elseif elitism
    EliteSourceParms(μ, first_half, second_half, η)
  elseif includeCenter
    CenterSelSourceParms(first_half, second_half)
  else
    SingleSourceParms(first_half, second_half)
  end
end

function SingleSourceParms(λ_1::Int, λ_2::Int)
    #order of vcat is very important
    o_orig = fill(:orig, λ_1)
    o_best = fill(:best, λ_2)
    SingleSourceParms(vcat(o_orig,o_best), SingleSource)
end
function AllSourceParms(μ::Int, λ_1::Int, λ_2::Int, η::Int)
  p_orig = fill((:orig,:parents), Int(floor(μ/η)))
  p_best = fill((:best,:parents), Int(floor(μ/η)))
  c_orig = fill((:orig,:center), 1)
  c_best = fill((:best,:center), 1)
  o_orig = fill((:orig,:offspring), λ_1)
  o_best = fill((:best,:offspring), λ_2)
  MultiSourceParms(vcat(p_orig,p_best,c_orig,c_best,o_orig,o_best), EliteCenterSource)
end

function EliteSourceParms(μ::Int, λ_1::Int, λ_2::Int, η::Int)
    p_orig = fill((:orig,:parents), Int(floor(μ/η)))
    p_best = fill((:best,:parents), Int(floor(μ/η)))
    o_orig = fill((:orig,:offspring), λ_1)
    o_best = fill((:best,:offspring), λ_2)
  MultiSourceParms(vcat(p_orig,p_best,o_orig,o_best), EliteSource)
end

function CenterSelSourceParms(λ_1::Int, λ_2::Int)
  c_orig = fill((:orig,:center), 1)
  c_best = fill((:best,:center), 1)
  o_orig = fill((:orig,:offspring), λ_1)
  o_best = fill((:best,:offspring), λ_2)
  MultiSourceParms(vcat(c_orig, c_best, o_orig, o_best), CenterSource)
end

# not used? - SelectionSource(ri::RunInfo, state::CMAES_State)                          = SelectionSource(sourceparms(ri), state) #sourceparms? shouldnt it be sourcevalues?
SelectionSource(sourceParms::SelectionSourceParms, state::CMAES_State)    = SelectionSource(sourceParms, sortorder(state))
#SelectionSource(sourceParms::SingleSourceParms, sortOrder::SortStructure) = SingleSource()
# have to change up SelectionSource to use η
function SelectionSource(sourceParms::SelectionSourceParms, sortOrder::SortStructure)
  sourceValues = sourcevalues(sourceParms)
  #println("SelectionSource_CMAES.jl::92 ->  sourcevalues before sort = $(sourceValues)")
  members = reshape(sourceValues, (1, length(sourceValues)))
  sPopn = SortedPopulation(members, sortOrder)
  source = vec(sPopn[:chr,:])
  #println("SelectionSource_CMAES.jl::96 -> source = $(source)")
  fit = sPopn[:fit, :]
  rank = tiedrank(fit)
  #println("SelectionSource_CMAES.jl::99 -> index = $(sortOrder.index)")
  selectiontype(sourceParms)(source, fit, rank)
end

combinedsize(s::SelectionSource) = length(s.fitness)

centerselected(s::SelectionSource)             = (any((x)->(x[2] == :center), s.source) && any((x)->(x[1] == :orig), s.source))

count(s::SelectionSource, origin::Symbol)      = count((x)-> (x[2] == origin), s.source)
ranks(s::SelectionSource, origin::Symbol)      = s.rank[s.source .== origin]
fitvalues(s::SelectionSource, origin::Symbol)  = s.fitness[s.source .== origin]

count(s::SelectionSource, source::Tuple)       = map((x::Symbol)->count(s, x), source)
ranks(s::SelectionSource, origins::Tuple)      = map((x::Symbol)->ranks(s, x), origins)
fitvalues(s::SelectionSource, source::Tuple)   = map((x::Symbol)->fitvalues(s, x), source)

parentfreq(sys::CMAES_System, pCount::Int)    = pCount/mu(sys)
offspringfreq(sys::CMAES_System, oCount::Int) = oCount/lambda(sys)

count(s::SelectionSource, species::Symbol, origin::Symbol) = count((x)-> (x[1] == species && x[2] == origin), s.source) # for any species of origin
count_species(s::SelectionSource, species::Symbol) = count((x)-> (x == species), s.source) # for single source
centerselected_(s::SelectionSource)            = any((x)->(x[2] == :center), s.source) && any((y)->(y[2] == :best),s.source)

function proportions(ecs::EliteCenterSource)
  totalCount = combinedsize(ecs)
  oCount = count(ecs, :offspring)						     # offspring count 	(Int)
  cCount_orig = centerselected(ecs) ? 1 : 0
  cCount_best = centerselected_(ecs) ? 1 : 0
  pCount = totalCount - oCount - cCount_orig - cCount_best         # parent count 		(Int)
  (pCount/totalCount, (cCount_orig + cCount_best)/totalCount, oCount/totalCount)
end

function proportions_species(ecs::EliteCenterSource)
    totalCount = combinedsize(ecs)
    oCount_orig = count(ecs, :orig, :offspring)
    oCount_best = count(ecs, :best, :offspring)
    cCount_orig = centerselected(ecs) ? 1 : 0
    cCount_best = centerselected_(ecs) ? 1 : 0
    pCount_orig = count(ecs, :orig, :parents)
    pCount_best = count(ecs, :best, :parents)
    (pCount_orig/totalCount, pCount_best/totalCount, cCount_orig/totalCount, cCount_best/totalCount, oCount_orig/totalCount, oCount_best/totalCount)
end

function proportions(es::EliteSource)
  totalCount = combinedsize(es)
  oCount = count(es, :offspring)
  pCount = totalCount - oCount
  (pCount/totalCount, NaN, oCount/totalCount)
end

function proportions_species(es::EliteSource)
    totalCount = combinedsize(es)
    oCount_orig = count(es, :orig, :offspring)
    oCount_best = count(es, :best, :offspring)
    pCount_orig = count(es, :orig, :parents)
    pCount_best = count(es, :best, :parents)
    (pCount_orig/totalCount, pCount_best/totalCount, NaN, NaN, oCount_orig/totalCount, oCount_best/totalCount)
end

function proportions(cs::CenterSource)
  totalCount = combinedsize(cs)
  oCount = count(cs, :offspring)
  cCount = centerselected(cs) ? 1 : 0
  (NaN, cCount/totalCount, oCount/totalCount)
end

function proportions_species(cs::CenterSource)
    totalCount = combinedsize(cs)
    oCount_orig = count(cs, :orig, :offspring)
    oCount_best = count(cs, :best, :offspring)
    cCount_orig = centerselected(cs) ? 1 : 0
    cCount_best = centerselected_(cs) ? 1 : 0
    (NaN, NaN, cCount_orig/totalCount, cCount_best/totalCount, oCount_orig/totalCount, oCount_best/totalCount)
end

function proportions(source::SingleSource)
    (NaN, NaN, NaN)
end

function proportions_species(source::SingleSource)
  totalCount = combinedsize(source)
  oCount_orig = count_species(source, :orig)
  oCount_best = count_species(source, :best)
  (NaN, NaN, NaN, NaN, oCount_orig/totalCount, oCount_best/totalCount)
end

avg_(sourceRank::Vector) = (length(sourceRank) > 0) ? mean(sourceRank) : NaN
max_(sourceRank::Vector) = (length(sourceRank) > 0) ? maximum(sourceRank) : NaN
min_(sourceRank::Vector) = (length(sourceRank) > 0) ? minimum(sourceRank) : NaN

rankmarginals(source::SelectionSource) = sourcemarginals(source, ranks)
fitmarginals(source::SelectionSource) = sourcemarginals(source, fitvalues)

function sourcemarginals(ecs::EliteCenterSource, fn::Function)
  parents, center, offspring = fn(ecs, (:parents, :center, :offspring))
  (max_(parents), avg_(parents), min_(parents), avg_(center),
   max_(offspring), avg_(offspring), min_(offspring))
end

function sourcemarginals(es::EliteSource, fn::Function)
  parents, offspring = fn(es, (:parents, :offspring))
  (max_(parents), avg_(parents), min_(parents), NaN,
   max_(offspring), avg_(offspring), min_(offspring))
end

function sourcemarginals(cs::CenterSource, fn::Function)
  center, offspring = fn(cs, (:center, :offspring))
  (NaN, NaN, NaN, avg_(center), max_(offspring), avg_(offspring), min_(offspring))
end

function sourcemarginals(source::SingleSource, fn::Function)
  (NaN, NaN, NaN, NaN, NaN, NaN, NaN)
end

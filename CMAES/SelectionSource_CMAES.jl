# struct SelectonSourceParms
#     sourceValues::Vector{Symbol}
#     selectionType::Type
# end

selectiontype(sParms::SelectionSourceParms) = sParms.selectionType
allsources(sParms::SelectionSourceParms) = sParms.sourcevalues          # depreciated
sourcevalues(sParms::SelectionSourceParms) = sParms.sourceValues

SelectionSourceParms(sys::CMAES_System) = SelectionSourceParms(elitism(sys), includecenter(sys), mu(sys), lambda(sys), η(sys))
AllSourceParms(sys::CMAES_System) = AllSourceParms(mu(sys), lambda(sys), η(sys))
EliteSourceParms(sys::CMAES_System) = EliteSourceParms(mu(sys), lambda(sys), η(sys))
CenterSelSourceParms(sys::CMAES_System) = CenterSelSourceParms(lambda(sys))

function SelectionSourceParms(elitism::Bool, includeCenter::Bool, μ::Int, λ::Int, η::Int)
  if elitism && includeCenter
    AllSourceParms(μ, λ, η)
  elseif elitism
    EliteSourceParms(μ, λ, η)
  elseif includeCenter
    CenterSelSourceParms(λ)
  else
    SingleSourceParms()
  end
end

function AllSourceParms(μ::Int, λ::Int, η::Int)
  p = fill(:parents, Int(floor(μ/η)))
  c = fill(:center, 1)
  o = fill(:offspring, λ)
  MultiSourceParms(vcat(p, c, o), EliteCenterSource)
end

function EliteSourceParms(μ::Int, λ::Int, η::Int)
  p = fill(:parents, Int(floor(μ/η)))
  o = fill(:offspring, λ)
  MultiSourceParms(vcat(p, o), EliteSource)
end

function CenterSelSourceParms(λ::Int)
  c = fill(:center, 1)
  o = fill(:offspring, λ)
  MultiSourceParms(vcat(c, o), CenterSource)
end


# abstract type SelectionSource end

# mutable struct EliteSource <: SelectionSource
#   source::Vector{Symbol}
#   fitness::Vector
#   rank::Vector{Float64}
# end

# mutable struct CenterSource <: SelectionSource
#   source::Vector{Symbol}
#   fitness::Vector
#   rank::Vector{Float64}
# end

# mutable struct EliteCenterSource <: SelectionSource
#   source::Vector{Symbol}
#   fitness::Vector
#   rank::Vector{Float64}
# end

# not used? - SelectionSource(ri::RunInfo, state::CMAES_State)                          = SelectionSource(sourceparms(ri), state) #sourceparms? shouldnt it be sourcevalues?
SelectionSource(sourceParms::SelectionSourceParms, state::CMAES_State)    = SelectionSource(sourceParms, sortorder(state))
SelectionSource(sourceParms::SingleSourceParms, sortOrder::SortStructure) = SingleSource()
# have to change up SelectionSource to use η
function SelectionSource(sourceParms::SelectionSourceParms, sortOrder::SortStructure)
  sourceValues = sourcevalues(sourceParms)
  #println("unsorted sourcevalues = $(sourceValues) of size $(length(sourceValues))")
  members = reshape(sourceValues, (1, length(sourceValues)))
  sPopn = SortedPopulation(members, sortOrder)
  source = vec(sPopn[:chr,:])
  #println("sorted sourcevalues = $(source) of size $(length(source))")
  fit = sPopn[:fit, :]
  rank = tiedrank(fit)
  selectiontype(sourceParms)(source, fit, rank)
end

combinedsize(s::SelectionSource) = length(s.fitness)

centerselected(s::SelectionSource)             = any((x)->(x == :center), s.source)
count(s::SelectionSource, origin::Symbol)      = count((x)-> (x == origin), s.source)
ranks(s::SelectionSource, origin::Symbol)      = s.rank[s.source .== origin]
fitvalues(s::SelectionSource, origin::Symbol)  = s.fitness[s.source .== origin]

count(s::SelectionSource, source::Tuple)       = map((x::Symbol)->count(s, x), source)
ranks(s::SelectionSource, origins::Tuple)      = map((x::Symbol)->ranks(s, x), origins)
fitvalues(s::SelectionSource, source::Tuple)   = map((x::Symbol)->fitvalues(s, x), source)

parentfreq(sys::CMAES_System, pCount::Int)    = pCount/mu(sys)
offspringfreq(sys::CMAES_System, oCount::Int) = oCount/lambda(sys)

function proportions(ecs::EliteCenterSource)
  totalCount = combinedsize(ecs)
  oCount = count(ecs, :offspring)						     # offspring count 	(Int)
  cCount   = centerselected(ecs) ? 1 : 0
  pCount = totalCount - oCount - cCount	         # parent count 		(Int)
  (pCount/totalCount, cCount/totalCount, oCount/totalCount)
end

function proportions(es::EliteSource)
  totalCount = combinedsize(es)
  oCount = count(es, :offspring)
  pCount = totalCount - oCount
  (pCount/totalCount, NaN, oCount/totalCount)
end

function proportions(cs::CenterSource)
  totalCount = combinedsize(cs)
  oCount = count(cs, :offspring)
  cCount = centerselected(cs) ? 1 : 0
  (NaN, cCount/totalCount, oCount/totalCount)
end

function proportions(source::SingleSource)
  (NaN, NaN, NaN)
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

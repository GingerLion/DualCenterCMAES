function sel☾μ▴λ☽(offspring::Population, nW::ShapedNoise, μ::Int)
  #sort and truncate one set of solutions whether offspring or offspring + parents
  (sMembers, sortOrder) = truncation(offspring, μ)
  sW = sort(nW, sortOrder)
  (sMembers, sW)
end

function sel☾μ▴λ✢1☽(offspring::Population, W::ShapedNoise, μ::Int, centerm::Member)
  center = RegularPopulation(offspring, centerm)
  allMembers = pcat(offspring, center)
  allW = pcat(W, CenterNoise(W))
  sel☾μ▴λ☽(allMembers, allW, μ)
end

function sel☾μ✢λ☽(c::CMAES_State, offspring::Population, nW::ShapedNoise, μ::Int; bug = false)
  segCount = segmentcount(c)
  ((segCount > 1) ? sel☾μ✢λ☽_(c, offspring, nW, μ, segCount, bug = bug)
  : (η(system(c)) >= 2) ? sel☾μ_η✢λ☽_(c, offspring, nW, μ, bug = bug)
  : sel☾μ✢λ☽_(c, offspring, nW, μ, bug = bug))
end

function sel☾μ✢λ✢1☽(c::CMAES_State, offspring::Population, nW::ShapedNoise, μ::Int, center::Member; bug = false)
  segCount = segmentcount(c)
  ((segCount > 1) ? sel☾μ✢λ✢1☽_(c, offspring, nW, μ, segCount, center, bug = bug) :
  (η(system(c)) >= 2) ? sel☾μ_η✢λ✢1☽_(c, offspring, nW, μ, center, bug = bug) :
  sel☾μ✢λ✢1☽_(c, offspring, nW, μ, center, bug = bug))
end

# We are inside an evaluation, so sOffspring has not been set yet this generation, so still holds the parents
# Between evaluations, pOffspring holds the parents
#(should consider changing so pOffspring set as parents at start of evolve instead of at end)
function parents✢noise(c::CMAES_State; bug = false)
    if bug
        (c.sOffspring_bug, ShapedNoise(c.sOffspring_bug, c.nModel_bug, bug = bug))
    else
        (c.sOffspring, ShapedNoise(c.sOffspring, c.nModel, bug = bug))
    end
end

function sel☾μ✢λ☽_(c::CMAES_State, offspring::Population, W::ShapedNoise, μ::Int; bug = false)
  (parents, pW) = parents✢noise(c, bug = bug)
  allMembers = pcat(parents, offspring)
  allW = pcat(pW, W)
  sel☾μ▴λ☽(allMembers, allW, μ)
end

function sel☾μ_η✢λ☽_(c::CMAES_State, offspring::Population, W::ShapedNoise, μ::Int; bug = false)
    #parents is sOffspring, pW is shaped noise of sOffspring
  (parents, pW) = parents✢noise(c, bug = bug)
  #update sortOrder of parents and return the sorted subset of the popn
  #(parents_new, sortOrder_new) = truncation(parents,Int(floor(μ/η(system(c)))))       #doesn't have to be μ/2 -> 2 will become a parameter in the system
  #population concatenate just the subset of parents and all offspring
  parents_new = RegularPopulation(members(parents)[:,1:Int(floor(μ/η(system(c))))],fitness(parents)[1:Int(floor(μ/η(system(c))))])
  allMembers = pcat(parents_new, offspring)
  pW_new = noisevalues(pW)[:chr,1:Int(floor(μ/η(system(c))))]
  allW = pcat(ShapedNoise(pW_new), W)
  sel☾μ▴λ☽(allMembers, allW, μ)
end

function sel☾μ✢λ✢1☽_(c::CMAES_State, offspring::Population, W::ShapedNoise, μ::Int, centerMember::Member; bug = false)
  (parents, pW) = parents✢noise(c, bug = bug)
  center = RegularPopulation(offspring, centerMember)
  allMembers = pcat(parents, center, offspring)
  allW = pcat(pW, CenterNoise(W), W)
  sel☾μ▴λ☽(allMembers, allW, μ)
end

function sel☾μ_η✢λ✢1☽_(c::CMAES_State, offspring::Population, W::ShapedNoise, μ::Int, centerMember::Member; bug = false)
  (parents, pW) = parents✢noise(c, bug = bug)
  #update sortOrder of parents and return the sorted subset of the popn
  #(parents_new, sortOrder_new) = truncation(parents,Int(floor(μ/η(system(c)))))       #doesn't have to be μ/2 -> 2 will become a parameter in the system
  center = RegularPopulation(offspring, centerMember)
  parents_new = RegularPopulation(members(parents)[:,1:Int(floor(μ/η(system(c))))],fitness(parents)[1:Int(floor(μ/η(system(c))))])
  allMembers = pcat(parents_new, center, offspring)
  pW_new = noisevalues(pW)[:chr,1:Int(floor(μ/η(system(c))))]
  allW = pcat(ShapedNoise(pW_new), CenterNoise(W), W)
  sel☾μ▴λ☽(allMembers, allW, μ)
end

function interleave_parms(c::CMAES_State)
  (c.sys.sParms.beginBinding, c.sys.sParms.loaded)
end

function sel☾μ✢λ☽_(c::CMAES_State, offspring::Population, W::ShapedNoise, μ::Int, segCount::Int)
  (parents, pW) = parents✢noise(c)
  (beginBinding, loaded) = interleave_parms(c)
  (allMembers, slotSize, pOrder, allOrder) = interleave(parents, offspring; beginBinding = beginBinding, loaded = loaded)
  cW = matchedinterleave(pW, W, slotSize, pOrder, allOrder)
  (sMembers, indexOrder) = truncation(allMembers, slotSize, μ, segCount)
  sW = matchedtruncation(cW, indexOrder)
  (sMembers, sW)
end

function sel☾μ✢λ✢1☽_(c::CMAES_State, nOffspring::Population, W::ShapedNoise, μ::Int, segCount::Int, centerm::Member; bug = false)
  error("μplusλ for add center and segmentation is not yet implemented")
end

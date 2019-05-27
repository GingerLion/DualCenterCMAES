function evolve_parms(state::CMAES_State)
  s = selection_parms(system(state))
  model = deepcopy(state.nModel)
  model_shadow = deepcopy(state.nModel_shadow)
  (state.gen, model, state.gen_shadow, model_shadow, s.N, s.λ, s.μ)
end

function update!(state::CMAES_State, model::CMAES_Model, model_shadow::CMAES_Model)
  state.pModel = state.nModel
  state.nModel = model
  state.pModel_shadow = state.nModel_shadow
  state.nModel_shadow = model_shadow
end

function update!(state::CMAES_State, nOffspring::Population, nOffspring_shadow::Population, sOffspring::Population, sOffspring_shadow::Population)
  state.pOffspring = state.sOffspring
  state.pOffspring_shadow = state.sOffspring_shadow
  state.nOffspring = nOffspring
  state.nOffspring_shadow = nOffspring_shadow
  state.sOffspring = sOffspring
  state.sOffspring_shadow = sOffspring_shadow
  #update shadow's center_ with best chromosome of sOffspring_shadow
  #center!_(state, bestchromosome(sOffspring_shadow))
end

function update!(state::CMAES_State, nE::Noise, nW::Noise, nW_shadow::Noise, sW::Noise, sW_shadow::Noise)
  state.pE = state.nE
  state.pW = state.sW
  state.pW_shadow = state.sW_shadow
  state.nE = nE
  state.nW = nW
  state.nW_shadow = nW_shadow
  state.sW = sW
  state.sW_shadow = sW_shadow
end

#updatebest!(state::CMAES_State) = (state.best = best(state.sOffspring))
function updatebest!(state::CMAES_State)
    if any(i -> isnan(i), members(state.sOffspring_shadow))
        state.best_shadow = ([NaN], NaN)
    else
        state.best = best(state.sOffspring)
        state.best_shadow = best(state.sOffspring_shadow)
        #println("bestfit = $(bestfitness(state.sOffspring)) \n bestfit_shadow = $(bestfitness(state.sOffspring_shadow))")
    end
end

function updategen!(state::CMAES_State)
  state.evalCount += evalsPerGen(system(state))
  if status_(state) != :found
      state.evalCount_shadow += evalsPerGen_(system(state))
      state.gen_shadow += 1
  else
      state.evalCount_shadow = state.evalCount_shadow
      state.gen_shadow = state.gen_shadow
      #println("not incrementing shadow system's evals.")
  end
  state.gen +=  1

end

function update!(state::CMAES_State, model::CMAES_Model, model_shadow::CMAES_Model,
                 nOffspring::Population, nOffspring_shadow::Population, sOffspring::Population, sOffspring_shadow::Population,
                 nE::Noise, nW::Noise, nW_shadow::Noise, sW::Noise, sW_shadow::Noise)
  # update state
  update!(state, model, model_shadow)
  update!(state, nOffspring, nOffspring_shadow, sOffspring, sOffspring_shadow)
  update!(state, nE, nW, nW_shadow, sW, sW_shadow)

  # update progress
  updategen!(state)
  updatebest!(state)
end

# -------------------------------
# Main Evolve rountine
# -------------------------------
# pFoo -> previous foo, nFoo -> new foo, sFoo -> selected foo

function evolvepopn!(state::CMAES_State, f::RealFitness)
  (gen, model, gen_shadow, model_shadow) = evolve_parms(state)

  #generate samples from center
  (nOffspring, nE, nW) = generatesamples(model)
  #(nOffspring_, nE_, nW_) = generatesamples_(model, state, deepcopy(nE))

  #generate shadow samples from center & and best chromosome
  #(nOffspring_shadow, nE, nW_shadow) = generatesamples(model_shadow, deepcopy(nE))
  (nOffspring_shadow, nOffspring_shadow_, nE_, nW_shadow, nW_shadow_) = generatesamples_(model_shadow, deepcopy(nE))   #generate samples based on best solution

  orig_λ = 0
  best_λ = 0
  (typeof(nW_shadow) <: Noise) ? orig_λ = popnsize(nW_shadow) : orig_λ = 0
  (typeof(nW_shadow_) <: Noise) ? best_λ = popnsize(nW_shadow_) : best_λ = 0

  #evaluate the popns separately
  evaluate!(nOffspring, f)
  if (typeof(nOffspring_shadow) <: Population)  evaluate!(nOffspring_shadow, f)  end
  if (typeof(nOffspring_shadow_) <: Population) evaluate!(nOffspring_shadow_, f) end

  #println("center_orig fitnesses :\n $(fitness(nOffspring_shadow))")
  #println("center_best fitnesses :\n $(fitness(nOffspring_shadow_))")
  #add shadow popn and shadow popn generated from the 2nd center and their fitness values
  if !(typeof(nOffspring_shadow) <: Population)
      dualcenterpopn = nOffspring_shadow_
      dualcenternoise = nW_shadow_
  elseif !(typeof(nOffspring_shadow_) <: Population)
      dualcenterpopn = nOffspring_shadow
      dualcenternoise = nW_shadow
  else
      dualcenterpopn = pcat(nOffspring_shadow, nOffspring_shadow_)
      #add shadow noise and shadow_ noise generated from the 2nd center
      dualcenternoise = nW_shadow + nW_shadow_
  end

  # select samples (and noise that generated them)
  (sOffspring, sW) = es_selection(state, model, f, nOffspring, nW)
  #(sOffspring_shadow, sW_shadow) = es_selection(state, model_shadow, f, nOffspring_shadow, nW_shadow, shadow = true)
  #println("dualcenterpopn members = $(members(dualcenterpopn))\n")
  #println("dualcenterpopn fitness = $(fitness(dualcenterpopn))\n")
  (sOffspring_shadow, sW_shadow) = es_selection(state, model_shadow, f, dualcenterpopn, dualcenternoise, shadow = true)
  #println("sortOrder = $(sOffspring_shadow.index)")
  #println("selected fitnesses = \n $(fitness(sOffspring_shadow))")

  flatten!(sW, weights(model))

  # fill array with :orig and :best symbols then sort them based on sortOrder and make them a population
  orig = fill(:orig, orig_λ)
  best = fill(:best, best_λ)
  all = vcat(orig, best)
  #println("Evolve_CMAES.jl::131 ->      source before sort = $(all)")
  memberz = reshape(all, (1, length(all)))
  sPopn = SortedPopulation(memberz, SortStructure(sOffspring_shadow))
  source = vec(sPopn[:chr,:])
  mems = deepcopy(members(sW_shadow))
  #println("μ = $(mu(state)), λ = $(lambda(state))")
  #println("source = $(source)")
  for i=1:length(source)
      if source[i] == :best
          mems[:,i] = (sOffspring_shadow[:chr, i] - center(model_shadow)) / sigma(model_shadow)
      end
  end
  #println("Evolve_CMAES.jl:142 -> source = $(source)")
  #println("center = $(center(model_shadow))")
  #println("sortOrder = $(sOffspring_shadow.sortOrder)")
  #println("Evolve_CMAES.jl::145 -> index = $(sOffspring_shadow.index)")
  #println("fitnesses = $(fitness(sOffspring_shadow))")
  #println("offspring = $(members(sOffspring_shadow))")
  #println("noise = $(mems)")
  #println("max euclidian distance at -> $(argmax(map(x -> norm(mems[:,x]), 1:length(source))))")
  #println("min euclidian distance at -> $(argmin(map(x -> norm(mems[:,x]), 1:length(source))))")
  sW_shadow = ShapedNoise(mems)

  flatten!(sW_shadow, weights(model_shadow))
  #flatten_with_σ!(sW_shadow, weights(model_shadow), model_shadow)

  nModel = update(model, sW, gen)
  nModel_shadow = update(model_shadow, sW_shadow, gen_shadow)

  update!(state, nModel, nModel_shadow, nOffspring, dualcenterpopn, sOffspring, sOffspring_shadow, nE, nW, dualcenternoise, sW, sW_shadow)

  generatesamples_prep!(state, nModel)
  if (any(i -> !isnan(i), nModel_shadow.C))
      generatesamples_prep!(state, nModel_shadow, shadow = true)
  end

  if status(state) != :found
      found!(state, f)
  end          # will overwrite status if found, even if the zero eigenvalue error was previously set
  if status_(state) != :found
      found!_(state,f)
  end

  if status(state) == :found
      status(state)
  elseif status_(state) == :found
      status_(state)             # returns the status
  end
end

function generatesamples_prep!(state::CMAES_State, model::CMAES_Model; shadow = false)
  eigendecomp!(state, model, shadow = shadow)
  sqrtC!(state, model, shadow = shadow)
  invsqrtC!(state, model, shadow = shadow)
end

function generatesamples(model::CMAES_Model)
  nE = SphericalNoise(chrlength(model), lambda(model))
  nW = ShapedNoise(nE, model)
  nOffspring = model + nW
  (nOffspring, nE, nW)
end
#shadowing generate sampels
#function generatesamples(model::CMAES_Model, nE)
#    nW = ShapedNoise(nE, model)
#    nOffspring = model + nW
#    (nOffspring, nE, nW)
#end
#shadowing generate samples off of best chromosome a.k.a center_
function generatesamples_(model::CMAES_Model, nE)
    (nW_orig, nW_best) = ShapedNoise(nE, model, dualcenter = true)
    !(typeof(nW_orig) <: Noise) ? nOffspring = NaN : nOffspring = model + nW_orig
    !(typeof(nW_best) <: Noise) ? nOffspring_best = NaN : nOffspring_best = RegularPopulation(nW_best, sigma(model) * members(nW_best) .+ center_(model))

    if !(typeof(nW_orig) <: Noise) && (typeof(nW_best) <: Noise)
        (NaN, nOffspring_best, nE, NaN, nW_best)
    elseif !(typeof(nW_best) <: Noise) && (typeof(nW_orig) <: Noise)
        (nOffspring, NaN , nE, nW_orig, NaN)
    else
        (nOffspring, nOffspring_best, nE, nW_orig, nW_best)
    end
end

function es_selection(state, model, f, nOffspring, nW; shadow = false)
  sys = system(state)
  μ = mu(sys)
  addparents = elitism(sys)
  addcentr = includecenter(sys)
  center = addcentr ? centermember(model, objfn(f)) : :nil
  center_ = addcentr ? centermember_(model,objfn(f)) : :nil
  #if any of the shadowged fitness goes from -Inf (to NaN) set status to found so both algorithms stop
  if any(i -> isnan(i), nOffspring[:chr,:])
      println("found NaNs in nOffspring")
      state.status = :found
      state.status_shadow = :found
      sOffspring = RegularPopulation(fill(NaN,(μ,size(nW[:chr,:],1))))
      nW.values = RegularPopulation(fill(NaN,(μ,size(nW[:chr,:],1))))
      return (sOffspring, nW)
  elseif shadow
      # when intializing, offspring == center & no parents or ν_population yet
      (initializing(state)    ? sel☾μ▴λ☽(nOffspring, nW, μ, shadow = true)                    :
       addparents && addcentr ? sel☾μ✢λ✢1☽(state, nOffspring, nW, μ, center, shadow = true) :
       addparents             ? sel☾μ✢λ☽(state, nOffspring, nW, μ, shadow = true)              :
       addcentr               ? sel☾μ▴λ✢1☽(nOffspring, nW, μ, center)
                              : sel☾μ▴λ☽(nOffspring, nW, μ, shadow = true))
  else
      # when intializing, offspring == center & no parents or ν_population yet
      (initializing(state)    ? sel☾μ▴λ☽(nOffspring, nW, μ)                    :
       addparents && addcentr ? sel☾μ✢λ✢1☽(state, nOffspring, nW, μ, center) :
       addparents             ? sel☾μ✢λ☽(state, nOffspring, nW, μ)              :
       addcentr               ? sel☾μ▴λ✢1☽(nOffspring, nW, μ, center)
                              : sel☾μ▴λ☽(nOffspring, nW, μ))

  end
end

function evolve_parms(state::CMAES_State)
  s = selection_parms(system(state))
  model = deepcopy(state.nModel)
  (state.gen, model, s.N, s.λ, s.μ)
end

function evolve_parms_(state::CMAES_State)
    s = selection_parms(system(state))
    model_shadow = deepcopy(state.nModel_shadow)
    (state.gen_shadow, model_shadow, s.N, s.λ, s.μ)
end

function update!(state::CMAES_State, model::CMAES_Model)
  state.pModel = state.nModel
  state.nModel = model
end

function update!_(state::CMAES_State, model::CMAES_Model)
    state.pModel_shadow = state.nModel_shadow
    state.nModel_shadow = model
end

function update!(state::CMAES_State, nOffspring::Population, sOffspring::Population)
  state.pOffspring = state.sOffspring
  state.nOffspring = nOffspring
  state.sOffspring = sOffspring
end

function update!_(state::CMAES_State, nOffspring::Population, sOffspring::Population)
  state.pOffspring_shadow = state.sOffspring_shadow
  state.nOffspring_shadow = nOffspring
  state.sOffspring_shadow = sOffspring
  #update shadow's center_ with best chromosome of sOffspring_shadow
  #center!_(state, bestchromosome(sOffspring_shadow))
end

function update!(state::CMAES_State, nE::Noise, nW::Noise, sW::Noise)
  state.pE = state.nE
  state.pW = state.sW
  state.nE = nE
  state.nW = nW
  state.sW = sW
end

function update!_(state::CMAES_State, nW::Noise, sW::Noise)
  state.pW_shadow = state.sW_shadow
  state.nW_shadow = nW
  state.sW_shadow = sW
end

#updatebest!(state::CMAES_State) = (state.best = best(state.sOffspring))
function updatebest!(state::CMAES_State)
    if any(i -> isnan(i), members(state.sOffspring_shadow))
        state.best_shadow = ([NaN], NaN)
    else
        state.best = best(state.sOffspring)
        state.best_shadow = best(state.sOffspring_shadow)
        if state.gen == 1
            state.best_overall = state.best
        end
        if state.gen_shadow == 1
            state.best_overall_ = state.best_shadow
        end
        if state.gen > 1
            if better_overall(state, state.best)
                state.best_overall = state.best
            end
        end
        if state.gen_shadow > 1
            if better_overall_(state, state.best_shadow)
                state.best_overall_ = state.best_shadow
            end
        end

        #println("bestfit = $(bestfitness(state.sOffspring)) \n bestfit_shadow = $(bestfitness(state.sOffspring_shadow))")
    end
end

function updategen!(state::CMAES_State)
  if !stopEvals(state) state.evalCount += evalsPerGen(system(state)); state.gen +=  1 end
  if !stopEvals_(state) state.evalCount_shadow += evalsPerGen_(system(state)); state.gen_shadow += 1 end
end

function update!(state::CMAES_State, model::CMAES_Model,
                 nOffspring::Population, sOffspring::Population,
                 nE::Noise, nW::Noise, sW::Noise)
  # update state
  update!(state, model)
  update!(state, nOffspring, sOffspring)
  update!(state, nE, nW, sW)

  # update progress
  updategen!(state)
  updatebest!(state)
end

function update!_(state::CMAES_State, model::CMAES_Model,
                 nOffspring::Population, sOffspring::Population, nW::Noise, sW::Noise)
  # update state
  update!_(state, model)
  update!_(state, nOffspring, sOffspring)
  update!_(state, nW, sW)

  # update progress
  updategen!(state)
  updatebest!(state)
end

# -------------------------------
# Main Evolve rountine
# -------------------------------
# pFoo -> previous foo, nFoo -> new foo, sFoo -> selected foo
function evolvepopn!(state::CMAES_State, f::RealFitness)
    (gen, model) = evolve_parms(state)

    (nOffspring, nE, nW) = generatesamples(model)

    evaluate!(nOffspring, f)

    # select samples (and noise that generated them)
    (sOffspring, sW) = es_selection(state, model, f, nOffspring, nW)

    flatten!(sW, weights(model))

    nModel = update(model, sW, gen)

    update!(state, nModel, nOffspring, sOffspring, nE, nW, sW)

    generatesamples_prep!(state, nModel)

    if status(state) != :found
        found!(state, f)
    end

    if status(state) == :found
        status(state)
    end
end
function evolvepopn!_(state::CMAES_State, f::RealFitness)
  (gen, model) = evolve_parms_(state) # do separate functions for this

  #generate shadow samples from center & and best chromosome
  #put this in other function (nE will be passed as function parameter or retrieved through the state)
  #if the other system isn't evolving then make it's own random noise
  if status(state) == :evolve
      (nOffspring, nOffspring_, nE, nW, nW_) = generatesamples_(model, state.nE)   #generate samples based on best solution
  else
      (nOffspring, nOffspring_, nE, nW, nW_) = generatesamples_(model, SphericalNoise(chrlength(state), lambda(system(state))))
  end
  orig_λ = 0
  best_λ = 0
  (typeof(nW) <: Noise) ? orig_λ = popnsize(nW) : orig_λ = 0
  (typeof(nW_) <: Noise) ? best_λ = popnsize(nW_) : best_λ = 0

  #evaluate the popns separately
  if (typeof(nOffspring) <: Population)  evaluate!(nOffspring, f)  end
  if (typeof(nOffspring_) <: Population) evaluate!(nOffspring_, f) end

  #println("center_orig fitnesses :\n $(fitness(nOffspring_shadow))")
  #println("center_best fitnesses :\n $(fitness(nOffspring_shadow_))")
  #add shadow popn and shadow popn generated from the 2nd center and their fitness values
  if !(typeof(nOffspring) <: Population)
      dualcenterpopn = nOffspring_
      dualcenternoise = nW_
  elseif !(typeof(nOffspring_) <: Population)
      dualcenterpopn = nOffspring
      dualcenternoise = nW
  else
      dualcenterpopn = pcat(nOffspring, nOffspring_)
      #add shadow noise and shadow_ noise generated from the 2nd center
      dualcenternoise = nW + nW_
  end

  #(sOffspring_shadow, sW_shadow) = es_selection(state, model_shadow, f, nOffspring_shadow, nW_shadow, shadow = true)
  #println("dualcenterpopn members = $(members(dualcenterpopn))\n")
  #println("dualcenterpopn fitness = $(fitness(dualcenterpopn))\n")
  (sOffspring, sW) = es_selection(state, model, f, dualcenterpopn, dualcenternoise, shadow = true)
  #println("sortOrder = $(sOffspring_shadow.index)")
  #println("selected fitnesses = \n $(fitness(sOffspring_shadow))")

  # fill array with :orig and :best symbols then sort them based on sortOrder and make them a population
  orig = fill(:orig, orig_λ)
  best = fill(:best, best_λ)
  all = vcat(orig, best)
  #println("Evolve_CMAES.jl::131 ->      source before sort = $(all)")
  memberz = reshape(all, (1, length(all)))
  sPopn = SortedPopulation(memberz, SortStructure(sOffspring))
  source = vec(sPopn[:chr,:])
  mems = deepcopy(members(sW))
  #println("μ = $(mu(state)), λ = $(lambda(state))")
  #println("source = $(source)")
  for i=1:length(source)
      if source[i] == :best
          mems[:,i] = (sOffspring[:chr, i] - center(model)) / sigma(model)
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
  sW = ShapedNoise(mems)

  flatten!(sW, weights(model))
  #flatten_with_σ!(sW_shadow, weights(model_shadow), model_shadow)

  nModel_ = update(model, sW, gen)

  update!_(state, nModel_, dualcenterpopn, sOffspring, dualcenternoise, sW)

  if (any(i -> !isnan(i), nModel_.C))
      generatesamples_prep!(state, nModel_, shadow = true)
  end

  if status_(state) != :found
      found!_(state,f)
  end

  if status_(state) == :found
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

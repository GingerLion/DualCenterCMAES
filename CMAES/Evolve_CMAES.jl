function evolve_parms(state::CMAES_State)
  s = selection_parms(system(state))
  model = deepcopy(state.nModel)
  model_bug = deepcopy(state.nModel_bug)
  (state.gen, model, state.gen_bug, model_bug, s.N, s.λ, s.μ)
end

function update!(state::CMAES_State, model::CMAES_Model, model_bug::CMAES_Model)
  state.pModel = state.nModel
  state.nModel = model
  state.pModel_bug = state.nModel_bug
  state.nModel_bug = model_bug
end

function update!(state::CMAES_State, nOffspring::Population, nOffspring_bug::Population, sOffspring::Population, sOffspring_bug::Population)
  state.pOffspring = state.sOffspring
  state.pOffspring_bug = state.sOffspring_bug
  state.nOffspring = nOffspring
  state.nOffspring_bug = nOffspring_bug
  state.sOffspring = sOffspring
  state.sOffspring_bug = sOffspring_bug
end

function update!(state::CMAES_State, nE::Noise, nW::Noise, nW_bug::Noise, sW::Noise, sW_bug::Noise)
  state.pE = state.nE
  state.pW = state.sW
  state.pW_bug = state.sW_bug
  state.nE = nE
  state.nW = nW
  state.nW_bug = nW_bug
  state.sW = sW
  state.sW_bug = sW_bug
end

#updatebest!(state::CMAES_State) = (state.best = best(state.sOffspring))
function updatebest!(state::CMAES_State)
    if any(i -> isnan(i), members(state.sOffspring_bug))
        state.best_bug = ([NaN], NaN)
    else
        state.best = best(state.sOffspring)
        state.best_bug = best(state.sOffspring_bug)
    end
end

function updategen!(state::CMAES_State)
  state.evalCount += evalsPerGen(system(state))
  state.evalCount_bug += evalsPerGen(system(state))
  state.gen +=  1
  state.gen_bug += 1
end

function update!(state::CMAES_State, model::CMAES_Model, model_bug::CMAES_Model,
                 nOffspring::Population, nOffspring_bug::Population, sOffspring::Population, sOffspring_bug::Population,
                 nE::Noise, nW::Noise, nW_bug::Noise, sW::Noise, sW_bug::Noise)
  # update state
  update!(state, model, model_bug)
  update!(state, nOffspring, nOffspring_bug, sOffspring, sOffspring_bug)
  update!(state, nE, nW, nW_bug, sW, sW_bug)

  # update progress
  updategen!(state)
  updatebest!(state)
end

# -------------------------------
# Main Evolve rountine
# -------------------------------
# pFoo -> previous foo, nFoo -> new foo, sFoo -> selected foo

function evolvepopn!(state::CMAES_State, f::RealFitness)
  (gen, model, gen_bug, model_bug) = evolve_parms(state)

  (nOffspring, nE, nW) = generatesamples(model)
  #use same nE
  (nOffspring_bug, nE, nW_bug) = generatesamples(model_bug, deepcopy(nE))
  #use it's own nE
  #(nOffspring_bug, nE, nW_bug) = generatesamples(model_bug)
  evaluate!(nOffspring, f)
  evaluate!(nOffspring_bug, f)

  # select samples (and noise that generated them)
  (sOffspring, sW) = es_selection(state, model, f, nOffspring, nW)
  (sOffspring_bug, sW_bug) = es_selection(state, model_bug, f, nOffspring_bug, nW_bug, bug = true)

  #w_1!(model_bug)
  flatten!(sW, weights(model))
  flatten!(sW_bug, weights(model_bug))

  nModel = update(model, sW, gen)
  nModel_bug = update(model_bug, sW_bug, gen_bug)

  update!(state, nModel, nModel_bug, nOffspring, nOffspring_bug, sOffspring, sOffspring_bug, nE, nW, nW_bug, sW, sW_bug)

  generatesamples_prep!(state, nModel)
  if (any(i -> !isnan(i), nModel_bug.C))
      generatesamples_prep!(state, nModel_bug)
  end

  found!(state, f)          # will overwrite status if found, even if the zero eigenvalue error was previously set
  found!_(state,f)
  if ((state.status != :found) & (state.status_bug == :found))
      println("shadow finished first")
      status_(state)
  else
      status(state)             # returns the status
  end
end

function generatesamples_prep!(state::CMAES_State, model::CMAES_Model)
  eigendecomp!(state, model)
  sqrtC!(state, model)
  invsqrtC!(state, model)
end

function generatesamples(model::CMAES_Model)
  nE = SphericalNoise(chrlength(model), lambda(model))
  nW = ShapedNoise(nE, model)
  nOffspring = model + nW
  (nOffspring, nE, nW)
end

function generatesamples(model::CMAES_Model, nE)
    nW = ShapedNoise(nE, model)
    nOffspring = model + nW
    (nOffspring, nE, nW)
end

function es_selection(state, model, f, nOffspring, nW; bug = false)
  sys = system(state)
  μ = mu(sys)
  addparents = elitism(sys)
  addcentr = includecenter(sys)
  center = addcentr ? centermember(model, objfn(f)) : :nil
  #if any of the bugged fitness goes from -Inf (to NaN) set status to found so both algorithms stop
  if any(i -> isnan(i), nOffspring[:chr,:])
      state.status = :found
      state.status_bug = :found
      sOffspring = RegularPopulation(fill(NaN,(μ,size(nW[:chr,:],1))))
      nW.values = RegularPopulation(fill(NaN,(μ,size(nW[:chr,:],1))))
      return (sOffspring, nW)
  elseif bug
      # when intializing, offspring == center & no parents or ν_population yet
      (initializing(state)    ? sel☾μ▴λ☽(nOffspring, nW, μ)                    :
       addparents && addcentr ? sel☾μ✢λ✢1☽(state, nOffspring, nW, μ, center, bug = true) :
       addparents             ? sel☾μ✢λ☽(state, nOffspring, nW, μ, bug = true)              :
       addcentr               ? sel☾μ▴λ✢1☽(nOffspring, nW, μ, center)
                              : sel☾μ▴λ☽(nOffspring, nW, μ))
  else
      # when intializing, offspring == center & no parents or ν_population yet
      (initializing(state)    ? sel☾μ▴λ☽(nOffspring, nW, μ)                    :
       addparents && addcentr ? sel☾μ✢λ✢1☽(state, nOffspring, nW, μ, center) :
       addparents             ? sel☾μ✢λ☽(state, nOffspring, nW, μ)              :
       addcentr               ? sel☾μ▴λ✢1☽(nOffspring, nW, μ, center)
                              : sel☾μ▴λ☽(nOffspring, nW, μ))

  end
end

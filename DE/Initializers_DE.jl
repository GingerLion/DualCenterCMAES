# struct PopnInitializer
#   initializers::Dict
#   function PopnInitializor(initMin::Array{Real}, initMax::Array{Real})
# end

function PopnInitializer(fit::DE_Fitness)
  pInit = PopnInitializer()
  pInit[:unif] = iunif_fn(fit.initMin, fit.initMax)
  pInit[:norm] = inorm_fn(fit.initMin, fit.initMax)
  pInit[:norm_de] = inormde_fn(fit.initMin, fit.initMax)
  pInit
end

function iunif_fn(iMin::Array{Float64,1}, iMax::Array{Float64,1})
  iDiff = iMax - iMin
  (de::DE_System) -> iMin .+ rand(de.chrLength, de.popnSize) .* iDiff
end

# normal distriubiton with the same mean and variance as the uniform distriubution
function inorm_fn(iMin::Array{Float64,1}, iMax::Array{Float64,1})			
  iDiff = iMax - iMin
  iAvg = (iMax + iMin) / 2
  (de::DE_System) -> iAvg .+ randn(de.chrLength, de.popnSize) / sqrt(12) .* iDiff
end

# clamped version of inorm from original julia code 
# which has about half the std dev of the unif distribution (sqrt(12)/6 = sqrt(1/3) = 0.577) 
function inormde_fn(iMin::Array{Float64,1}, iMax::Array{Float64,1})
  iDiff = iMax - iMin
  iAvg = (iMax + iMin) / 2
  (de::DE_System) -> iAvg .+ clamp.(randn(de.chrLength, de.popnSize) / 6, -0.5 , 0.5) .* iDiff
end

getindex(p::PopnInitializer, key::Symbol) = p.initializers[key]
setindex!(p::PopnInitializer, value::Function, key::Symbol) = (p.initializers[key] = value)


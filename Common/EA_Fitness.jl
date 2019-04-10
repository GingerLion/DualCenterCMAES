# struct UnBoundedFitness{T} <: EA_Fitness
#   objFn::Function
#   dimension::Int64
#   direction::Symbol     # either :min or :max respectively when minimizing or maximizing the function
#   optimalValue::T    	  # if known, else set to Inf or -Inf depending on direction
# 	optimum::Vector       # if known, else set to [NaN, ..., NaN]
#   epsilon::T            # difference from goal that is considered a success
# end

# struct BoundedFitness{T} <: EA_Fitness
#   objFn::Function
#   dimension::Int64
#   direction::Symbol     # either :min or :max respectively when minimizing or maximizing the function
#   optimalValue::T    	  # if know else set to Inf or -Inf depending on direction
# 	optimum::Vector       # if known, else set to [NaN, ..., NaN]
#   epsilon::T         	  # difference from goal that is considered a success
#   initMax::Array{T}  	  # upper bound on the solution vector
#   initMin::Array{T}  	  # lower bound on the solution vector
# end

function Fitness(objFn::Function, dimensions::Integer, direction::Symbol, epsilon;
					optimalValue = (direction == :min) ? -Inf : Inf,
					optimum = fill(dimensions, NaN))
  UnBoundedFitness(objFn, dimensions, direction, optimalValue, optimum, epsilon)
end

function Fitness(objFn::Function, dimension::Integer, direction::Symbol, optimalValue::T, optimum::Vector,
						epsilon::T, initMax::Vector, initMin::Vector) where {T}
#					   optimalValue = (direction == :min) ? -Inf : Inf,
#					   optimum = fill(dimensions, NaN))
  BoundedFitness(objFn, dimension, direction, optimalValue, optimum, epsilon, initMax, initMin)
end

function Fitness(objFn::Function, dimension::Integer, direction::Symbol, optimalValue::T, optimum::Vector,
						epsilon::T, initMax::Float64, initMin::Float64) where {T}
					   # optimalValue = (direction == :min) ? -Inf : Inf, 
					   # optimum = fill(dimensions, NaN))
  BoundedFitness(objFn, dimension, direction, optimalValue, optimum, epsilon, fill(initMax, dimension), fill(initMin, dimension))
end


returntype(f::UnBoundedFitness{T}) where {T} = T
returntype(f::BoundedFitness{T})   where {T} = T

objfn(f::Fitness)		 = f.objFn
dimensions(f::Fitness) 	 = f.dimension
direction(f::Fitness) 	 = f.direction
optimalvalue(f::Fitness) = f.optimalValue
optimum(f::Fitness) 	 = f.optimum

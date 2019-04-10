# mutable struct BestFitHistory
#   history::Vector{Float64}
#   ptr::Int
#   windowSize::Int
#   direction::Symbol
#   function BestFitHistory(windowSize::Int, direction::Symbol) 
#   	bfh= new()
#   	bfh.windowSize = windowSize
#   	bfh.history = fill(Inf, windowSize)
# 	bfh.ptr = 1
# 	bfh.direction = direction
# 	bfh
#   end
# end

length(bfh::BestFitHistory) = bfh.windowSize

function setindex!(bfh::BestFitHistory, value::Float64)
	bfh.history[bfh.ptr] = value
	bfh.ptr = bfh.ptr % length(bfh) + 1
end

maximizing(bfh::BestFitHistory) = bfh.direction == :max
minimizing(bfh::BestFitHistory) = bfh.direction == :min
maximum(bfh::BestFitHistory) = maximum(bfh.history)
minimum(bfh::BestFitHistory) = minimum(bfh.history)
best(bfh::BestFitHistory) = (maximizing(bfh) ? maximum(bfh) : minimum(bfh))
worst(bfh::BestFitHistory) = (minimizing(bfh) ? maximum(bfh) : minimum(bfh))
history(bfh::BestFitHistory) = copy(bfh.history)
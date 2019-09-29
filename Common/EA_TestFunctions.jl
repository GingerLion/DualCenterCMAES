#---------------------------------------------
# Linear Transformation
#---------------------------------------------

pairrotate(α, i, j) = givens(cos(α),sin(α),i,j)[1]

function Rotation(d::Int, α)
	diag = fill(1.0, d)
	rotatn = (Diagonal(diag) + zeros(d,d))
	for i = d:-1:2
		rotatn = pairrotate(α, i-1, i) * rotatn
	end
	rotatn
end

Rotation˚(d::Int, α˚) = Rotation(d, deg2rad(α˚))

function condrot˚(d::Int, α˚)
	M = Rotation(d, deg2rad(α˚))
	(eigval, eigvec) = eig(M)
#	maximum(eigval) / minimum(eigval)
end

function lineartransform(fn, d, α˚, Δ)
	M = Rotation˚(d, α˚)
	(x) -> fn(M * (x - Δ))
end

function lineartransform(fn, d, α˚)
	M = Rotation˚(d, α˚)
	(x) -> fn(M * x)
end

#---------------------------------------------
# Test Function Definitions
#---------------------------------------------
# VALLEY SHAPED FUNCTIONS (rosenbrock, dixon_price)
# Rosenbrock :
rosenbrock(x) = sum(map((i) -> ((1.0 - x[i])^2 + 100.0 * (x[i+1] - x[i]^2)^2), 1:(length(x)-1)))
rosenbrock() = typeof(rosenbrock([0.0, 0.0]))
# Dixon-Price :
function dixon_price(x)
	  term1 = (x[1] - 1)^2
	  term2 = sum(map((i)->(i*(2x[i]^2 - x[i-1])^2),2:length(x)))
	  term1 + term2
end
# Rastrigin :

rastrigin(x) =  10.0 * length(x) + sum(x.^2 - 10.0cos.(2π * x))

rastrigin_(x) = rastrigin(x - ones(length(x)))

rastrigin() = typeof(rastrigin([0.0 ,0.0]))



# PLATE SHAPED FUNCTIONS
# zakharov :
function zakharov(x)
       term1 = sum(x.^2)
       term2 = sum(map((i)->(0.5*i*x[i]),1:length(x)))^2
       term3 = sum(map((i)->(0.5*i*x[i]),1:length(x)))^4
       term1 + term2 + term3
end

# Swechfel

schwefel(x) =  (418.982887 * length(x)) - sum(x .* sin.(sqrt.(abs.(x))))
schwefel() = typeof(rastrigin([0.0 ,0.0]))

function schwefel(x)
  d = length(x)
  total = sum(x .* sin.(sqrt.(abs.(x))))
  418.982887 * d - total
end
# Sphere

sphere(x) = sum(x .* x)
sphere() = typeof(sphere([0.0, 0.0]))

# Elliptical
function elliptical(x)
  sum_ = 0
  for i = 1:length(x)
  	for j = 1:i
  		sum_ += x[j]^2
  	end
  end
  sum_
end

# Ackley
function ackley(x)
	d = length(x)
    a = 20.0
    b = 0.2
    c = 2π

	sum1 = sum(x .^ 2)
	sum2 = sum(cos.(c * x))

	term1 = -a * exp.(-b * sqrt.(sum1/d))
	term2 = -exp.(sum2 / d)

	term1 + term2 + a + exp(1)
end

# Levy
function levy(x)
	d = length(x)
	w = 1.0 .+ (x .- 1.0) / 4

	term1 = (sin(π * w[1]))^2
	term3 = (w[d]-1.0)^2 * (1.0 + sin(2π*w[d])^2)

	wi = w[1:(d-1)]
	total = sum((wi .- 1.0).^2 .* (1.0 .+ 10.0 * sin.(π * wi .+ 1.0).^2))

	term1 + total + term3
end

# griewank
function griewank(x)
  i = collect(1:length(x))
  total = sum(x .^ 2 / 4000)
  product = prod(cos.(x ./ sqrt.(i)))

  total - product + 1
end

griewank() = typeof(griewank([0.0, 0.0]))

#---------------------------------------------
# Test Function Generators
#---------------------------------------------

function generatetests(n; ε = 1.0e-10)
	allZeros = fill(0.0, n)
	allOnes = fill(1.0, n)
	EA_Test = Dict()
	EA_Test[:ackley] = UnBoundedFitness{Float64}(ackley, n, :min, 0.0, allZeros, ε)
	EA_Test[:levy] = UnBoundedFitness{Float64}(levy, n, :min, 0.0, allOnes, ε)
	EA_Test[:griewank] = UnBoundedFitness{Float64}(griewank, n, :min, 0.0, allZeros, ε)
	EA_Test[:rastrigin] = UnBoundedFitness{Float64}(rastrigin, n, :min, 0.0, allZeros, ε)
	EA_Test[:rastrigin_] = UnBoundedFitness{Float64}(rastrigin_, n, :min, 0.0, allZeros, ε)
	EA_Test[:rosenbrock] = UnBoundedFitness{Float64}(rosenbrock, n, :min, 0.0, allOnes, ε)
	EA_Test[:schwefel] = UnBoundedFitness{Float64}(schwefel, n, :min, 0.0, allZeros, ε)
	EA_Test[:sphere] = UnBoundedFitness{Float64}(sphere, n, :min, 0.0, allZeros, ε)
	EA_Test[:elliptical] = UnBoundedFitness{Float64}(elliptical, n, :min, 0.0, allZeros, ε)
	EA_Test
end

function boundedtests(n; ε = 1.0e-10)
	allZeros = fill(0.0, n)
	allOnes = fill(1.0, n)
	EA_Test = Dict()
	EA_Test[:ackley] = BoundedFitness{Float64}(ackley, n, :min, 0.0, allZeros, ε, fill(32.768, n), fill(-32.768, n))
	EA_Test[:levy] = BoundedFitness{Float64}(levy, n, :min, 0.0, allOnes, ε, fill(10.0, n), fill(-10.0, n))
	EA_Test[:griewank] = BoundedFitness{Float64}(griewank, n, :min, 0.0, allZeros, ε, fill(600.0, n), fill(-600.0, n))
	EA_Test[:rastrigin] = BoundedFitness{Float64}(rastrigin, n, :min, 0.0, allZeros, ε, fill(5.12, n), fill(-5.12, n))
	EA_Test[:rosenbrock] = BoundedFitness{Float64}(rosenbrock, n, :min, 0.0, allOnes, ε, fill(sqrt(3), n), fill(-sqrt(3), n))
	EA_Test[:schwefel] = BoundedFitness{Float64}(schwefel, n, :min, 0.0, allZeros, ε, fill(512.0, n), fill(-512.0, n))
	EA_Test[:sphere] = BoundedFitness{Float64}(sphere, n, :min, 0.0, allZeros, ε, fill(5.12, n), fill(-5.12, n))
	EA_Test[:elliptical] = BoundedFitness{Float64}(elliptical, n, :min, 0.0, allZeros, ε, fill(65.536, n), fill(-65.536, n))
	EA_Test
end

function generatetests(n, α˚; ε = 1.0e-10)
	allZeros = fill(0.0, n)
	allOnes = fill(1.0, n)
	EA_Test = Dict()
	EA_Test[:ackley] = UnBoundedFitness{Float64}(lineartransform(ackley, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test[:levy] = UnBoundedFitness{Float64}(lineartransform(levy, n, α˚), n, :min, 0.0, allOnes, ε)
	EA_Test[:griewank] = UnBoundedFitness{Float64}(lineartransform(griewank, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test[:rastrigin] = UnBoundedFitness{Float64}(lineartransform(rastrigin, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test[:rosenbrock] = UnBoundedFitness{Float64}(lineartransform(rosenbrock, n, α˚), n, :min, 0.0, allOnes, ε)
	EA_Test[:dixon_price] = UnBoundedFitness{Float64}(lineartransform(dixon_price, n, α˚), n, :min, 0.0,  map((i)->(2^-((2^i - 2)/2^i)), 1:n), ε)
	EA_Test[:zakharov] = UnBoundedFitness{Float64}(lineartransform(zakharov, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test[:rastrigin_] = UnBoundedFitness{Float64}(rastrigin_, n, :min, 0.0, allZeros, ε)
	EA_Test[:schwefel] = UnBoundedFitness{Float64}(lineartransform(schwefel, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test[:sphere] = UnBoundedFitness{Float64}(lineartransform(sphere, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test[:elliptical] = UnBoundedFitness{Float64}(lineartransform(elliptical, n, α˚), n, :min, 0.0, allZeros, ε)
	EA_Test
end

function boundedtests(n, α˚; ε = 1.0e-10)
	allZeros = fill(0.0, n)
	allOnes = fill(1.0, n)
	EA_Test = Dict()
	EA_Test[:ackley] = BoundedFitness{Float64}(lineartransform(ackley, n, α˚), n, :min, 0.0, allZeros, ε, fill(32.768, n), fill(-32.768, n))
	EA_Test[:levy] = BoundedFitness{Float64}(lineartransform(levy, n, α˚), n, :min, 0.0, allOnes, ε, fill(10.0, n), fill(-10.0, n))
	EA_Test[:griewank] = BoundedFitness{Float64}(lineartransform(griewank, n, α˚), n, :min, 0.0, allZeros, ε, fill(600.0, n), fill(-600.0, n))
	EA_Test[:rastrigin] = BoundedFitness{Float64}(lineartransform(rastrigin, n, α˚), n, :min, 0.0, allZeros, ε, fill(5.12, n), fill(-5.12, n))
	EA_Test[:rosenbrock] = BoundedFitness{Float64}(lineartransform(rosenbrock, n, α˚), n, :min, 0.0, allOnes, ε, fill(sqrt(3), n), fill(-sqrt(3), n))
	EA_Test[:schwefel] = BoundedFitness{Float64}(lineartransform(schwefel, n, α˚), n, :min, 0.0, allZeros, ε, fill(512.0, n), fill(-512.0, n))
	EA_Test[:sphere] = BoundedFitness{Float64}(lineartransform(sphere, n, α˚), n, :min, 0.0, allZeros, ε, fill(5.12, n), fill(-5.12, n))
	EA_Test[:elliptical] = BoundedFitness{Float64}(lineartransform(elliptical, n, α˚), n, :min, 0.0, allZeros, ε, fill(65.536, n), fill(-65.536, n))
	EA_Test
end

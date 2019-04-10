p1 = RegularPopulation([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0], direction = :min)
p2 = RegularPopulation(mcat(p1,p1), direction = :min)
p3 = p1 + p1
p1[1]
p2[:chr,1]
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
#rosenbrock() = typeof(rosenbrock([0,0]))
f = RealFitness(rosenbrock, 8, :min, 0.0, fill(0.0,8), 1e-10)

evaluate(p1, f)
evaluate!(p1, f)
p4 = pcat(p1,p1)
p4[:chr,2]
p4[:fit,2]
mbr = [1 2; 3 4; 5 6; 7 8;  9 10; 11 12; 13 14; 15 16]
simplefitfn(x) = abs(mean(x) - 8)
f = RealFitness(simplefitfn, 8, :min, 0.0, fill(0.0,4), 1e-10)

p5 = RegularPopulation(mbr)
evaluate!(p5, f)
s5 = SortedPopulation(p5)
s5[1]
s5[2]
s5[:fit,1]
s5[:fit,2]
s5[:fit,3]
s5[:chr,1]
s5[:chr,2]
s5[:chr,3]
srtOdr = truncate!(s5,2)
p6 = RegularPopulation([1 2 3 4 5])
s6 = SortedPopulation(p6, srtOdr)
s6[3]

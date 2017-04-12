include("2d-ising-monte-carlo.jl")

T = 15
L = 5
N1 = 100
N2 = 10000

#Compile
iterate_over_temperatures(Array(linspace(0.6, 0.8, 5)), 5, 5, 10)

#Time
tic()
iterate_over_temperatures(Array(linspace(0.6, 0.8, T)), L, N1, N2)
x = toq();

println("Spin flips per second is : ", L^2*(N1+N2)*T/1e6/x)

include("funcs.jl")
include("metropolis.jl")
include("wolf.jl")
include("wolf-spins.jl")
#include("plotting.jl")

T = 0.5
N = 100
N_eq = 10
L = 5

S = SpinArray(L, random_config(L))
neighbors(4,5, S.L)

cluster = falses(L,L)


i,j = 2,3
createcluster!(S, cluster, i, j, T)

S_old = copy(S.configuration)

flipcluster!(S, cluster)


S.configuration

cluster
S.configuration


S.configuration - S_old

#E, phi = zeros(N), zeros(N)
#flips = rand(1:L, 2*L^2, N+N_eq, 2)
#randoms = rand(Float32, 2*L^2, N+N_eq)


#metropolis_step!(S, Exp, flips[1,1,1],flips[1,1,2], randoms[1,1], T)

#metropolis_sweep!(S, Exp, flips[:,1,:], randoms[:,1], T)


E, M, Cv, X = iterate_over_temperatures(linspace(0.5, 4.0, 50), 6, 10, 500)

Ew, Mw, Cvw, Xw = iterate_over_temperatures_wolff(linspace(0.5, 4.0, 20), 6, 10, 500)

display(E)
display(Ew)

#IsingPlots.allplots(linspace(0.5, 4.0, 50), L, E, M, Cv, X)

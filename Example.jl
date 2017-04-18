# EXAMPLE
include("Metropolis2DIsing.jl")
include("IsingPlots.jl")

const L = 10
const N_eq = 100
const N_mc = 5000

T = Array(linspace(0.1, 5.0, 50))
E, M, Cv, X = Metropolis2DIsing.iterate_over_temperatures(T, L, N_eq, N_mc)
IsingPlots.allplots(T, L, E, M, Cv, X)

#Genatic Ising Funciton
type SpinArray
  L::Int #lattice size
  configuration::Array{Int8, 2} #lattice of +/-1 Ising spins
end

flip!(A::Array{Int8, 2}, i::Int, j::Int) = A[i, j] *= -1
random_config(L::Int) = rand(-1:2:1, L, L)
magnetization(A::Array{Int8, 2}) =  abs(sum(A))
periodic_index(a::Int, L::Int) = ( a == 0 ? L : (a == L + 1 ? 1 : a))

function energy_single_flip(spinarray::SpinArray, i::Int, j::Int)::Int #energy of configuration
    L = spinarray.L
    M = spinarray.configuration
    @inbounds return -2*M[i,j]*( M[i, periodic_index(j+1, L)] + M[i, periodic_index(j-1, L)] + M[periodic_index(i+1, L), j] + M[periodic_index(i-1, L), j] )
end

function energy(spinarray::SpinArray)::Int #energy of configuration
    E = 0
    L = spinarray.L
    M = spinarray.configuration
    @inbounds for j = 1:L, i = 1:L  #Sum the energy of the four nearest-neighbour spins
        E -= M[i,j]*(M[i, periodic_index(j+1, L)] + M[periodic_index(i+1, L), j] )
    end
    return E
end

function thermo_quantities(T::Float64, L::Int64, N_eq::Int64, N::Int64)::Tuple{Float64,Float64,Float64,Float64}

    S = SpinArray(L, random_config(L))
    Exp = [exp(-4.0/T), exp(-8.0/T)]
    E, phi = zeros(N), zeros(N)
    flips = rand(1:L, 2*L^2, 2, N+N_eq)
    randoms = rand(Float32, 2*L^2, N+N_eq)

    @inbounds for i = 1:N_eq #Run a few times to equilibriate
        metropolis_sweep!(S, Exp, flips[:,:,i], randoms[:,i], T)
    end

     @inbounds for i = 1:N #Runs which will be included in the thermodynamic averages
        j = i + N_eq
        metropolis_sweep!(S, Exp, flips[:,:,j], randoms[:,j], T)
        E[i] = energy(S)
        phi[i] = magnetization(S.configuration)
    end

    E_avg = sum(E)/(N*L^2) #average energy at T
    M = sum(phi)/(N*L^2) #average Mag at T
    Cv = (dot(E,E)/N - sum(E)^2/N^2)/(T^2)
    X = (dot(phi, phi)/N - sum(phi)^2/N^2)/(T)

    return E_avg, M, Cv, X
end

function iterate_over_temperatures(Ts, L::Int, N_eq::Int, N::Int)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    num_of_Ts = length(Ts)
    E, M, Cv, X = Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts)

    @inbounds for j = 1:length(Ts)
         E[j], M[j], Cv[j], X[j] = thermo_quantities(Ts[j], L, N_eq, N)
    end

    return E, M, Cv, X
end

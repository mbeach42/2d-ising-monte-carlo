
@everywhere function random_config(L::Int64)::Array{Int64,2} #generate a random (L,L) matrix of +/-1
    return 2*rand(0:1, L, L) - 1
end

@everywhere function shiftm(i::Int64, L::Int64)
  return (i+L-2)%L + 1
end
@everywhere function shiftp(i::Int64, L::Int64)
  return (i)%L + 1
end

@everywhere function magnetization(M::Array{Int64,2})::Int64
    return sum(M)
end

@everywhere function energy_singleflip(M::Array{Int64,2}, i::Int64, j::Int64, L::Int64)::Int64 #energy of configuration
    @inbounds temp = 2*M[i,j]*(M[i, shiftp(j,L)] + M[i, shiftm(j,L)] + M[shiftp(i,L), j] + M[shiftm(i,L), j])
    return temp
end

@everywhere function energy(M::Array{Int64,2}, L::Int64)::Int64 #energy of configuration
    E = 0
    for i = 1:L
        for j = 1:L #Sum the energy of the four nearest-neighbour spins
            @inbounds E += M[i,j]*(M[i, shiftp(j,L)] + M[shiftp(i,L), j])
        end
    end
    return div(E, 8)
end

@everywhere function MC_step!(config::Array{Int64,2},L::Int64, T::Float64, flipx::Int64, flipy::Int64,r::Float64 ) #Monte Carlo step
    @inbounds config[flipx, flipy] *= -1 #randomly flip one spin
    deltaE = 1.*energy_singleflip(config, flipx, flipy, L) #compute energy of spin flip
    if deltaE > 0 && r > exp(-1.0/T*deltaE) #if spin flip is un-favourable, revert back
        @inbounds config[flipx, flipy] *= -1 #reverse the spin flip
    end
end

@everywhere function MC_sweep!(config::Array{Int64,2}, T::Float64, L::Int64, flipx::Int64, flipy::Int64, r::Float64) #Sweep over L^2 MC steps
    for i = 1:2*L^2
        MC_step!(config, L,  T, flipx, flipy, r)
    end
end

@everywhere function thermo_quantities(T::Float64, L::Int64, N_eq::Int64, N_steps::Int64)::Tuple{Float64,Float64,Float64,Float64}

    config = random_config(L)
    E, phi = zeros(N_steps), zeros(N_steps)
    flipsx, flipsy, rs = rand(1:L, N_steps+N_eq), rand(1:L, N_steps+N_eq), rand(N_steps+N_eq)

    for i = 1:N_eq #Run a few times to equilibriate
        MC_sweep!(config, T, L, flipsx[i], flipsy[i], rs[i])
    end

    for i = 1:N_steps #Runs which will be included in the thermodynamic averages
        MC_sweep!(config, T, L, flipsx[i], flipsy[i], rs[i])
        E[i] = energy(config, L)
        phi[i] = magnetization(config)
    end

    E_avg = sum(E)/(N_steps*L^2) #average energy at T
    M = sum(phi)/(N_steps*L^2) #average Mag at T
    Cv = (dot(E,E)/N_steps - sum(E)^2/N_steps^2)/(T^2)
    X = (dot(phi, phi)/N_steps - sum(phi)^2/N_steps^2)/(T)
    return E_avg, M, Cv, X
end

@everywhere function iterate_over_temperatures(Ts::Array{Float64,1}, L::Int64,
    N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1}, Array{Float64,1},Array{Float64,1},Array{Float64,1}}

    F = length(Ts)
    E, M, Cv, X = zeros(F), zeros(F), zeros(F), zeros(F)

    for j = 1:length(Ts)
        E[j], M[j], Cv[j], X[j] = thermo_quantities(Ts[j], L, N_eq, N_steps)
    end
    return E, M, Cv, X
end

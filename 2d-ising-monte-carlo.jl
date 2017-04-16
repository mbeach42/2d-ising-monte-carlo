# A Metropolis Algorithm for the 2D Ising model
# Currently limited by rand() and the time it takes to do energy_singleflip

#using RandomNumbers.Xorshifts
#r = Xoroshiro128Plus();
#set_zero_subnormals(true)

@fastmath @everywhere function random_config(L::Int64)::Array{Int64,2} #generate a random (L,L) matrix of +/-1
    return rand(-1:2:1, L, L)
end

@everywhere function periodic_index(a::Int64, L::Int64)::Int64 #This impliments periodic boudnary conditions (PBC's)
    if a == 0
        return L
    elseif a == L + 1
        return 1
    else
      return a
  end
end

@fastmath @everywhere function magnetization(M::Array{Int64,2})::Int64
    @inbounds return abs(sum(M))
end

@fastmath @everywhere function energy_singleflip(M::Array{Int64,2},
                                                 i::Int64, j::Int64, L::Int64)::Int64 #energy of configuration
    @inbounds return -2*M[i,j]*( M[i, periodic_index(j+1, L)] + M[i, periodic_index(j-1, L)] + M[periodic_index(i+1, L), j] + M[periodic_index(i-1, L), j] )
end

@fastmath @everywhere function energy(M::Array{Int64,2}, L::Int64)::Int64 #energy of configuration
    E = 0
    @inbounds for j = 1:L, i = 1:L  #Sum the energy of the four nearest-neighbour spins
        E -= M[i,j]*(M[i, periodic_index(j+1, L)] + M[periodic_index(i+1, L), j] )
    end
    return E
end

@fastmath @everywhere function MC_step!(config::Array{Int64,2}, L::Int64, T::Float64,
                                        pre_exp::Array{Float64, 1},
                                        flipx::Int64, flipy::Int64, rand_num::Float64) #Monte Carlo step
    config[flipx, flipy] *= -1 #randomly flip one spin
    deltaE = 1.*energy_singleflip(config, flipx, flipy, L) #compute energy of spin flip
    if deltaE == 4
      rand_num > pre_exp[1] ? config[flipx, flipy] *= -1 : nothing
    elseif deltaE == 8
      rand_num > pre_exp[2] ? config[flipx, flipy] *= -1 : nothing
    end
    return nothing
end

@fastmath @everywhere function MC_sweep!(config::Array{Int64,2}, T::Float64, L::Int64,
                                        pre_exp::Array{Float64, 1},
                                        flipx::Array{Int64, 1}, flipy::Array{Int64, 1},
                                        rand_num::Array{Float64, 1}) #Sweep over L^2 MC steps
    @inbounds for i = 1:2*L^2
         MC_step!(config, L,  T, pre_exp, flipx[i], flipy[i], rand_num[i])
    end
    return nothing
end

@fastmath @everywhere function thermo_quantities(T::Float64, L::Int64,
                                                 N_eq::Int64, N_steps::Int64)::Tuple{Float64,Float64,Float64,Float64}
    config = random_config(L)
    E, phi = zeros(N_steps), zeros(N_steps)
    flipsx, flipsy, rs = rand(1:L, 2*L^2, N_steps+N_eq), rand(1:L, 2*L^2, N_steps+N_eq), rand(2*L^2, N_steps+N_eq)
    pre_exp = [exp(-4.0/T), exp(-8.0/T)]

    @inbounds for i = 1:N_eq #Run a few times to equilibriate
        MC_sweep!(config, T, L, pre_exp, flipsx[:, i], flipsy[:, i], rs[:, i])
    end
     @inbounds for i = 1:N_steps #Runs which will be included in the thermodynamic averages
        MC_sweep!(config, T, L, pre_exp, flipsx[:, N_eq+i], flipsy[:, N_eq+i], rs[:, N_eq+i])
        E[i] = energy(config, L)
        phi[i] = magnetization(config)
    end
    E_avg = sum(E)/(N_steps*L^2) #average energy at T
    M = sum(phi)/(N_steps*L^2) #average Mag at T
    Cv = (dot(E,E)/N_steps - sum(E)^2/N_steps^2)/(T^2)
    X = (dot(phi, phi)/N_steps - sum(phi)^2/N_steps^2)/(T)
    return E_avg, M, Cv, X
end

@fastmath @everywhere function iterate_over_temperatures(Ts::Array{Float64,1}, L::Int64,
                                                         N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    F = length(Ts)
    E, M, Cv, X = zeros(F), zeros(F), zeros(F), zeros(F)
    @inbounds for j = 1:length(Ts)
         E[j], M[j], Cv[j], X[j] = thermo_quantities(Ts[j], L, N_eq, N_steps)
    end
    return E, M, Cv, X
end

@fastmath @everywhere function iterate_over_temperatures_parallel(Ts::Array{Float64,1},
                                                                  L::Int64,
                                                                  N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    F = length(Ts)
    E, M, Cv, X = SharedArray{Float64}(F),  SharedArray{Float64}(F), SharedArray{Float64}(F), SharedArray{Float64}(F)
    @inbounds @sync @parallel for j = 1:length(Ts)
         E[j], M[j], Cv[j], X[j] = thermo_quantities(Ts[j], L, N_eq, N_steps)
    end
    return E, M, Cv, X
end

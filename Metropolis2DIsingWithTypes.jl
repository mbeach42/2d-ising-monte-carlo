# A Metropolis Algorithm for the 2D Ising model
# Currently limited by rand() and the time it takes to do energy_singleflip

module Metropolis2DIsingWithTypes

export iterate_over_temperatures, iterate_over_temperatures_parallel


immutable Config
    temp::Float64
    size::Int
    spins::Array{Int64, 2}
end

function random_config(L::Int64)::Array{Int64,2} #generate a random (L,L) matrix of +/-1
    return rand(-1:2:1, L, L)
end

function periodic_index(a::Int64, L::Int64)::Int64 #This impliments periodic boudnary conditions (PBC's)
    if a == 0
        return L
    elseif a == L + 1
        return 1
    else
      return a
  end
end

function magnetization(config::Config)::Int64
    return abs(sum(config.spins))
end

function energy_singleflip(config::Config, i::Int64, j::Int64)::Int64  #energy of configuration
    M, L = config.spins, config.size
    @inbounds return -2*M[i,j]*( M[i, periodic_index(j+1, L)] + M[i, periodic_index(j-1, L)] + M[periodic_index(i+1, L), j] + M[periodic_index(i-1, L), j] )
end

function energy(config::Config)::Int64 #energy of configuration
    M, L = config.spins, config.size
    E = 0.0
    @inbounds for j = 1:L, i = 1:L  #Sum the energy of the four nearest-neighbour spins
        E -= M[i,j]*(M[i, periodic_index(j+1, L)] + M[periodic_index(i+1, L), j] )
    end
    return E
end

function mcstep!(config::Config, pre_exp::Array{Float64, 1}, flipx::Int64, flipy::Int64, rand_num::Float64) #Monte Carlo step
    M, L, T = config.spins, config.size, config.temp
    M[flipx, flipy] *= -1 #randomly flip one spin
    deltaE = 1.*energy_singleflip(config, flipx, flipy) #compute energy of spin flip
    if deltaE == 4
      rand_num > pre_exp[1] ? M[flipx, flipy] *= -1 : nothing
    elseif deltaE == 8
      rand_num > pre_exp[2] ? M[flipx, flipy] *= -1 : nothing
    end
    return nothing
end

function mcsweep!(config::Config, pre_exp::Array{Float64, 1},
                  flipx::Array{Int64, 1}, flipy::Array{Int64, 1}, rand_num::Array{Float64, 1}) #Sweep over L^2 MC steps
    @inbounds for i = 1:2*config.size^2
         mcstep!(config, pre_exp, flipx[i], flipy[i], rand_num[i])
    end
    return nothing
end

function thermo_quantities(T::Float64, L::Int64, N_eq::Int64, N_steps::Int64)::Tuple{Float64,Float64,Float64,Float64}
    config = Config(T, L, random_config(L))
    E, phi = zeros(N_steps), zeros(N_steps)
    flipsx, flipsy, rs = rand(1:L, 2*L^2, N_steps+N_eq), rand(1:L, 2*L^2, N_steps+N_eq), rand(2*L^2, N_steps+N_eq)
    pre_exp = [exp(-4.0/T), exp(-8.0/T)]

    @inbounds for i = 1:N_eq #Run a few times to equilibriate
        mcsweep!(config, pre_exp, flipsx[:, i], flipsy[:, i], rs[:, i])
    end
     @inbounds for i = 1:N_steps #Runs which will be included in the thermodynamic averages
        mcsweep!(config, pre_exp, flipsx[:, N_eq+i], flipsy[:, N_eq+i], rs[:, N_eq+i])
        E[i] = energy(config)
        phi[i] = magnetization(config)
    end
    E_avg = sum(E)/(N_steps*L^2) #average energy at T
    M = sum(phi)/(N_steps*L^2) #average Mag at T
    Cv = (dot(E,E)/N_steps - sum(E)^2/N_steps^2)/(T^2)
    X = (dot(phi, phi)/N_steps - sum(phi)^2/N_steps^2)/(T)
    return E_avg, M, Cv, X
end

function iterate_over_temperatures(Ts::Array{Float64,1}, L::Int64, N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    F = length(Ts)
    E, M, Cv, X = zeros(F), zeros(F), zeros(F), zeros(F)
    @inbounds for j = 1:length(Ts)
         E[j], M[j], Cv[j], X[j] = thermo_quantities(Ts[j], L, N_eq, N_steps)
    end
    return E, M, Cv, X
end

function iterate_over_temperatures_parallel(Ts::Array{Float64,1}, L::Int64, N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    F = length(Ts)
    E, M, Cv, X = SharedArray{Float64}(F),  SharedArray{Float64}(F), SharedArray{Float64}(F), SharedArray{Float64}(F)
    @inbounds @sync @parallel for j = 1:length(Ts)
         E[j], M[j], Cv[j], X[j] = thermo_quantities(Ts[j], L, N_eq, N_steps)
    end
    return E, M, Cv, X
end


end

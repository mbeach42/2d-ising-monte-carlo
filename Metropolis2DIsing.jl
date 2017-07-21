using PyPlot
using HDF5
using NPZ
using ProfileView

# Types
type SpinArray
  L::Int                            # lattice size
  configuration::Array{Float64, 2}  # lattice of q-spin clock states
end

# Short Functions
random_config(L::Int, q::Int) = rand(1:q, L, L)/float(q)
magnetization(A::Array{Float64, 2}) =  abs( sum(sinpi.(A))^2 + sum(cospi.(A))^2 )/size(A)[1]^2.
periodic_index(a::Int, L::Int) = ( a == 0 ? L : (a == L + 1 ? 1 : a))

# Functions 
function neighbors(i, j, L)
    up = [periodic_index(i+1,L), j]
    down = [periodic_index(i-1,L), j]
    left = [i, periodic_index(j-1,L)]
    right = [i, periodic_index(j+1,L)]
    return up, down, left, right
end

function energy(spinarray::SpinArray)
    E = 0.0
    L = spinarray.L
    M = spinarray.configuration
    for j = 1:L, i = 1:L                 # Sum the energy of the four nearest-neighbour spins
        for (n,m) in neighbors(i, j, L)
            E -= cospi(M[i,j]-M[n,m])
        end
    end
    return E/4.0
end

function stiffness_x(spinarray::SpinArray, T::Float64)
    term1, term2 = 0.0, 0.0
    L = spinarray.L
    M = spinarray.configuration
    for j = 1:L, i = 1:L
        x = M[i,j] - M[i, periodic_index(j+1,L)]
        term1 += cospi(x)
        term2 += sinpi(x)
   end
   stiffness =  term1 - 1./T*term2^2
   return stiffness *= 1./(L^2)
end

function createcluster!(S::SpinArray, cluster::BitArray{2}, i::Int, j::Int, T::Float64)
    L = S.L
    cluster[i,j] = true
    for (n,m) in neighbors(i, j, L)
        if S.configuration[n,m] == S.configuration[i,j] && cluster[n,m] == false 
            rand() < 1. - exp(-2./T) ? createcluster!(S, cluster, n, m, T) : nothing   
        end
    end
end

function flipcluster!(S::SpinArray, cluster::BitArray{2}, q::Int)
    L = S.L
    rand_spin = rand(1:q)/float(q)
    for i in find(cluster)
        S.configuration[i] += 1/float(q)
    end
     S.configuration =  S.configuration.%1
end

function wolff!(S::SpinArray, T::Float64, q::Int, maxit::Int = 10)
    L = S.L
    for iter = 1:maxit
        cluster = falses(L,L)
        x, y = rand(1:L), rand(1:L)
        createcluster!(S, cluster, x, y, T)
        flipcluster!(S, cluster, q)
        #maximum(S.configuration) == minimum(S.configuration) ? break : nothing # break the loop if it's already found the solution
    end
end

function thermo_quantities_wolff(T::Float64, L = 10, N_eq = 10, N = 10, num_cluster = 10, q = 3, T_num = 1, configs= zeros(100,100))#::Tuple{Float64,Float64,Float64,Float64}
    S = SpinArray(L, random_config(L, q))
    E, phi, stiff = zeros(N), zeros(N), zeros(N)
    for k = 1:N_eq  # Warm-up
        wolff!(S, T, q, num_cluster)
     end
    for j = 1:N 
        wolff!(S, T, q, num_cluster)
        configs[(T_num-1)*N + j  , :] = reshape(S.configuration, L*L)
        E[j] = energy(S)
        phi[j] = magnetization(S.configuration)
        stiff[j] = stiffness_x(S, T)
    end
    E_avg = sum(E)/(N*L^2)          
    M = sum(phi)/(N*L^2)                       
    Cv = (dot(E,E)/N - sum(E)^2/N^2)/(T^2)
    X = (dot(phi, phi)/N - sum(phi)^2/N^2)/(T)
    Stiff = sum(stiff)/N
    return E_avg, M, Cv, X, Stiff
end

function iterate_over_temperatures_wolff(Ts; L = 10, N_eq = 10, N = 10, num_cluster = 10, q = 3)#::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    num_of_Ts = length(Ts)
    E, M, Cv, X, Stiff = Array{Float64}(num_of_Ts), Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts)
    
    configs = zeros(num_of_Ts*N, L*L)
    ispath("L=$L-q=$q-data/") == false ? mkdir("L=$L-q=$q-data/") : nothing # Create directory if it doesnt exist
    
    @inbounds @fastmath for j = 1:num_of_Ts
         E[j], M[j], Cv[j], X[j], Stiff[j] = thermo_quantities_wolff(Ts[j], L, N_eq, N, num_cluster, q, j, configs)
    end
    
    npzwrite("L=$L-q=$q-data/configs.npy", configs)

    return E, M, Cv, X, Stiff
end

# Plotting functions
function allplots(name, q, T_range, L, E, M, Cv, X)
    fig = figure("all_plot",figsize=(10,4))
    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".", label = "L = $L_val")
    legend()
    
    subplot(222)
    xlabel("Temperature")
    ylabel("Magnetization")
    plot(T_range, abs.(M), ".", color = "red")
    
    subplot(223)
    xlabel("Temperature")
    ylabel("Cv")
    plot(T_range, Cv, "." , color = "orange")
    
    subplot(224)
    xlabel("Temperature")
    ylabel("Sus")
    plot(T_range, X, ".", color = "green");
    
    savefig("L=$L-q=$q-data/"*name*".png")
end


#Parameters
T_range = linspace(0.2, 3.0, 40)
L_val = 28
N_eq = 100
N = 1000
num_cluster = 2000

#Q = 2
@time E, M, Cv, X, Stiff = iterate_over_temperatures_wolff(T_range; L = L_val, N_eq = N_eq, N = N, q = 2, num_cluster = num_cluster);
allplots("4plots", 2, T_range, L_val, E, M, Cv, X)

#Q = 6
@time E, M, Cv, X, Stiff = iterate_over_temperatures_wolff(T_range; L = L_val, N_eq = N_eq, N = N, q = 6, num_cluster = num_cluster);
allplots("4plots", 6, T_range, L_val, E, M, Cv, X)


#Q = 8
@time E, M, Cv, X, Stiff = iterate_over_temperatures_wolff(T_range; L = L_val, N_eq = N_eq, N = N, q = 8, num_cluster = num_cluster);
allplots("4plots", 8, T_range, L_val, E, M, Cv, X)
function thermo_quantities_wolff(T::Float64, L::Int64, N_eq::Int64, N::Int64)::Tuple{Float64,Float64,Float64,Float64}

    S = SpinArray(L, random_config(L))
    Exp = [exp(-4.0/T), exp(-8.0/T)]
    E, phi = zeros(N), zeros(N)
    flips = rand(1:L, 2*L^2, N+N_eq, 2)
    randoms = rand(Float32, 2*L^2, N+N_eq)

    @inbounds for i = 1:N_eq #Run a few times to equilibriate
        wolff!(S, T, 200)
    end

     @inbounds for i = 1:N #Runs which will be included in the thermodynamic averages
        j = i + N_eq
        wolff!(S, T, 200)
        E[i] = energy(S)
        phi[i] = magnetization(S.configuration)
    end

    E_avg = sum(E)/(N*L^2) #average energy at T
    M = sum(phi)/(N*L^2) #average Mag at T
    Cv = (dot(E,E)/N - sum(E)^2/N^2)/(T^2)
    X = (dot(phi, phi)/N - sum(phi)^2/N^2)/(T)

    return E_avg, M, Cv, X
end

function iterate_over_temperatures_wolff(Ts, L::Int, N_eq::Int, N::Int)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},Array{Float64,1}}
    num_of_Ts = length(Ts)
    E, M, Cv, X = Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts),Array{Float64}(num_of_Ts)

    @inbounds for j = 1:length(Ts)
         E[j], M[j], Cv[j], X[j] = thermo_quantities_wolff(Ts[j], L, N_eq, N)
    end

    return E, M, Cv, X
end

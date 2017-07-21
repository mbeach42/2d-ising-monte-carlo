#Metropolis Algorithm functions

function metropolis_step!(spinarray::SpinArray, precomputed_exp::Array{Float64, 1}, i::Int, j::Int64, r::Float32, T::Float64)

    flip!(spinarray.configuration, i, j) #randomly flip one spin

    ΔE = energy_single_flip(spinarray, i, j) #compute energy of a spin flip

    if ΔE == 4 #If the energy is positive, only accept the new spin if rand() < exp(-beta/T)
      r > precomputed_exp[1] ? flip!(spinarray.configuration, i, j) : nothing
    elseif ΔE == 8
      r > precomputed_exp[2] ? flip!(spinarray.configuration, i, j) : nothing
    end
end

function metropolis_sweep!(spinarray::SpinArray, precomputed_exp::Array{Float64, 1}, fliparray::Array{Int64, 2}, random_numbers::Array{Float32, 1}, T::Float64)

    @inbounds for i = 1:2*spinarray.L^2
         metropolis_step!(spinarray, precomputed_exp, fliparray[i, 1], fliparray[i, 2], random_numbers[i], T)
    end
end

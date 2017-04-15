include("2d-ising-monte-carlo.jl")
#include("MC_medium.jl")

#Compile
iterate_over_temperatures(Array(linspace(0.6, 0.8, 5)), 5, 5, 10)
iterate_over_temperatures_parallel(Array(linspace(0.6, 4, 5)), 5, 10, 50)

function timeit(T, L, N)
  time = zeros(10)
  for i = 1:10
    #Time
    tic()
    iterate_over_temperatures(Array(linspace(0.1, 4.0, T)), L, 10, N)
    time[i] = 2*L^2*(10+N)*T/1e6/toq();

  end

  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

function timeit_parallel(T, L, N)
  time = zeros(10)
  for i = 1:10
    #Time
    tic()
    iterate_over_temperatures_parallel(Array(linspace(0.1, 4.0, T)), L, 10, N)
    time[i] = L^2*(10+N)*T/1e6/toq();

  end

  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

include("Metropolis2DIsing.jl")

#Compile once
Metropolis2DIsing.iterate_over_temperatures(Array(linspace(0.6, 0.8, 5)), 5, 5, 10)

function timeit(T, L, N)
  time = zeros(10)
  for i = 1:10
    tic()
    Metropolis2DIsing.iterate_over_temperatures(Array(linspace(0.1, 4.0, T)), L, 10, N)
    time[i] = 2*L^2*(10+N)*T/1e6/toq();

  end
  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

timeit(20, 10, 1000)

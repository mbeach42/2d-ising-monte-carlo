include("2d-ising-monte-carlo.jl")

#Compile
iterate_over_temperatures(Array(linspace(0.6, 0.8, 5)), 5, 5, 10)


function timeit(T, L, N)
  time = zeros(100)
  for i = 1:100
    #Time
    tic()
    iterate_over_temperatures(Array(linspace(0.6, 0.8, T)), L, 10, N)
    time[i] = L^2*(10+N)*T/1e6/toq();

  end

  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

function timeit_parallel(T, L, N)
  time = zeros(100)
  for i = 1:100
    #Time
    tic()
    iterate_over_temperatures_parallel(Array(linspace(0.6, 0.8, T)), L, 10, N)
    time[i] = L^2*(10+N)*T/1e6/toq();

  end

  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

#Wolf Cluster Algorithm

function neighbors(i, j, L)
  up = [periodic_index(i+1,L), j]
  down = [periodic_index(i-1,L), j]
  left = [i, periodic_index(j-1,L)]
  right = [i, periodic_index(j+1,L)]
  return up, down, left, right
end

function createcluster!(spinarray::SpinArray, cluster::BitArray{2}, i::Int, j::Int, T::Float64)
    cluster[i,j] = true
    for (n,m) in neighbors(i, j, L)
        if S.configuration[n,m] == S.configuration[i,j] && cluster[n,m] == false #not sure if the order matters
           if rand() < 1. - exp(-2./T)
             createcluster!(spinarray, cluster, n, m, T)
           end
        end
    end
end

function flipcluster!(spinarray::SpinArray, cluster::BitArray{2})
    for i in find(cluster)
        spinarray.configuration[i] *= -1
    end
end


function wolff!(spinarray::SpinArray, T::Float64, maxit::Int = 300)

    for iter in 1:maxit
        cluster = falses(L,L)
        x, y = rand(1:L), rand(1:L)
        createcluster!(spinarray, cluster, x, y, T)
    end
    flipcluster!(spinarray, cluster)
end

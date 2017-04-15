function laplacian_bad(n)
  x = rand(n)
  @inbounds for i = 1:n
    a = x[i]
    b = x[i].^6  + sin(9)
  end
end

function laplacian_good(n)

  for i = 1:n
    @inbounds a = x[i]
    @inbounds b = x[i].^6  + sin(9)
  end
end

using Polyhedra
using TropicalFrechetMeans

function trop_normalize(x)
    return x .- first(x)
end

sample = [[-3,0,0], [0,-6,0], [0,0,-12]]

num_FM = tropical_frechet_mean(sample) |> trop_normalize
@show num_FM

P = tropical_frechet_set(sample; tol=1e-3)
@show P
@show vrep(P)
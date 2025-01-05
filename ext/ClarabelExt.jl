module ClarabelExt
using TropicalFrechetMeans, Polyhedra
using Clarabel, JuMP

TropicalFrechetMeans.polyhedral_frechet_model(sample, alphas; power=2) =
    polyhedral_frechet_model(Clarabel.Optimizer, sample, alphas; power=power)
TropicalFrechetMeans.polyhedral_frechet_mean(sample, alphas; power=2) = 
    polyhedral_frechet_mean(Clarabel.Optimizer, sample, alphas; power=power)

function TropicalFrechetMeans.polyhedral_frechet_set(lib::Lib, sample, alphas; power=2, tol=1e-3) where {
    Lib<:Polyhedra.Library
}
    H = polyhedral_frechet_set(Clarabel.Optimizer, lib, sample, alphas; power=power, tol=tol) |> hrep
    redH = removehredundancy(H, Clarabel.Optimizer)
    return polyhedron(redH, lib)
end

function TropicalFrechetMeans.tropical_frechet_set(lib::Lib, sample::Vector{Vector{T}}; power=2, tol=1e-3) where {
    Lib<:Polyhedra.Library, T<:Real
}
    H = tropical_frechet_set(Clarabel.Optimizer, lib, sample; power=power, tol=tol) |> hrep
    redH = removehredundancy(H, Clarabel.Optimizer)
    return polyhedron(redH, lib)
end

end # module ClarabelExt
module ClarabelExt # Should be same name as the file (just like a normal package)

using TropicalFrechetMeans, Clarabel, Polyhedra

TropicalFrechetMeans.polyhedral_frechet_model(sample, alphas; power=2, tol=1e-3) = 
    polyhedral_frechet_model(Clarabel.Optimizer, sample, alphas; power=power, tol=tol)

TropicalFrechetMeans.polyhedral_frechet_mean(sample, alphas; power=2, tol=1e-3) = 
    polyhedral_frechet_mean(Clarabel.Optimizer, sample, alphas; power=power, tol=tol)

TropicalFrechetMeans.tropical_frechet_mean(sample; power=2, tol=1e-3) = 
    tropical_frechet_mean(Clarabel.Optimizer, sample; power=power, tol=tol) 

function TropicalFrechetMeans.polyhedral_frechet_set(
    lib::Lib, 
    sample, 
    alphas; 
    power=2, tol=1e-3
) where {
    Lib<:Polyhedra.Library
}
    H = polyhedral_frechet_set(Clarabel.Optimizer, lib, sample, alphas; power=power, tol=tol) |> hrep
    
    redH = removehredundancy(H, () -> Clarabel.Optimizer(verbose=false))
    return polyhedron(redH, lib)
end

function TropicalFrechetMeans.tropical_frechet_set(lib::Lib, sample::Vector{Vector{T}}; power=2, tol=1e-3) where {
    Lib<:Polyhedra.Library, T<:Real
}
    return tropical_frechet_set(Clarabel.Optimizer, lib, sample; power=power, tol=tol)
end

end # module
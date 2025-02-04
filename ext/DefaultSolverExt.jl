module DefaultSolverExt # Should be same name as the file (just like a normal package)

using TropicalFrechetMeans, CDDLib, Clarabel

TropicalFrechetMeans.polyhedral_frechet_set(sample, alphas; power=2, tol=1e-3) = 
    polyhedral_frechet_set(CDDLib.Library(:exact), sample, alphas; power=power, tol=tol)

function TropicalFrechetMeans.tropical_frechet_set(sample; power=2, tol=1e-3)
    alphas = tropical_ball_facets(length(sample |> first))
    return polyhedral_frechet_set(sample, alphas; power=power, tol=tol)
end

end # module
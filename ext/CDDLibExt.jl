module CDDLibExt
using TropicalFrechetMeans, Polyhedra
using CDDLib

TropicalFrechetMeans.polyhedral_frechet_set(sample, alphas; power=2, tol=1e-3) = 
    polyhedral_frechet_set(CDDLib.Library(:exact), sample, alphas; power=power, tol=tol)

TropicalFrechetMeans.tropical_frechet_set(sample; power=2, tol=1e-3) = 
    tropical_frechet_set(CDDLib.Library(:exact), sample, alphas; power=power, tol=tol)

end # module ClarabelExt
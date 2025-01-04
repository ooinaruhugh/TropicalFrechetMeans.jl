using MathOptInterface
using LinearAlgebra
using JuMP, Polyhedra

# poly_frechet_set(sample, alphas; power=2, rep=:polyhedron, tol=1e-3) = poly_frechet_set(CDDLib.Library, sample, alphas; power=power, rep=rep, tol=tol)
"""
Find the set of polyhedral Fréchet means of a given sample.
The rows of `alphas` are the facet normals scaled to α⋅x = 1.
`power` gives the exponent of the distance before taking the sum.
`rep` is either "vrep" or "hrep" -- returning either vertices or halfspaces.
"""
function polyhedral_frechet_set(::Type{Opt}, ::Type{Lib}, sample, alphas; power=2, rep=:polyhedron, tol=1e-3) where {Opt<:MathOptInterface.AbstractOptimizer,Lib<:Polyhedra.Library}
    num_facets = size(alphas)[1]
    
    # Compute one Fréchet mean
    FM = polyhedral_frechet_mean(Opt, sample, alphas, power=power)
    
    @debug "Frechet mean found: $(FM)"
    
    # Rationalise all coordinates
    rat_alphas = rationalize.(alphas, tol=tol)
    rat_sample = [rationalize.(pt, tol=tol) for pt in sample]
    rat_mean = rationalize.(FM, tol=tol)
    
    distances = [polyhedral_distance(rat_mean, pt, alphas) for pt in rat_sample]
    
    Amat = vcat(rat_alphas, -rat_alphas)
    bval1 = Rational{Int64}[]
    
    evals = rat_alphas * hcat([rat_mean - pt for pt in rat_sample]...)
    for k in 1:length(rat_sample)
        evals[:,k] .-= distances[k]
    end
    
    progress = 0
    total = num_facets
    
    for k = 1:num_facets
        
        greatest_nonpos = argmax([x > 0 ? -Inf : x for x in evals[k,:]])
        push!(bval1, dot(rat_alphas[k,:], rat_sample[greatest_nonpos]) + 
            distances[greatest_nonpos])
        
        progress += 1
        @debug "Removing redundant half-spaces: $(round(progress/total * 100, digits=3))%   \r"    
    end
    
    @debug "\nFinding defining facets..."
    
    poly = polyhedron(hrep(rat_alphas, bval1), Lib(:exact))
    # poly = polyhedron(hrep(Amat, [bval1...,bval1...]), Lib(:exact))
    removehredundancy!(poly)
    
    if rep == :hrep
        @debug "Finding facets..."
        return hrep(poly)
    elseif rep == :vrep
        @debug "Finding vertices..."
        return vrep(poly)
    else
        @debug "Defaulting to polyhedron."
        return poly
    end
end


"""
    tropical_frechet_set(sample; power=2, rep=:polyhedron, tol=1e-3)

Find one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals scaled to α⋅x = 1.
`power` gives the exponent of the distance before taking the sum.
"""
function tropical_frechet_set(::Type{Opt}, ::Type{Lib}, sample; power=2, rep=:polyhedron, tol=1e-3) where {Opt<:MathOptInterface.AbstractOptimizer,Lib<:Polyhedra.Library}
    dim = length(sample[1])
    alphas = trop_ball_facets(dim)
    return polyhedral_frechet_set(Opt, Lib, sample, alphas; power=power, rep=rep, tol=tol)
end
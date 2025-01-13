using LinearAlgebra
using MathOptInterface
using JuMP
using Polyhedra
using Clarabel, CDDLib

"""
Find the set of polyhedral Fréchet means of a given sample.
The rows of `alphas` are the facet normals scaled to α⋅x = 1.
`power` gives the exponent of the distance before taking the sum.
`rep` is either "vrep" or "hrep" -- returning either vertices or halfspaces.
"""
function polyhedral_frechet_set(::Type{Opt}, lib::Lib, sample, alphas; power=2, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer,
    Lib<:Polyhedra.Library
}
    FM = if tol > 0 
        sample = map(sample) do pt; rationalize.(pt; tol=tol) end
        pt = polyhedral_frechet_mean(Opt, sample, alphas; power=power)
        rationalize.(pt; tol=tol)
    else 
        polyhedral_frechet_mean(Opt, sample, alphas; power=power)
    end
    @debug "Frechet mean found: $(FM)"

    B = map(sample) do pt 
        d = polyhedral_distance(FM, pt, alphas)
        polyhedral_ball(lib, pt, alphas, d)
    end
    
    return intersect(B...)
end
function polyhedral_frechet_set(::Type{Opt}, sample::Vector{Vector{T}}, alphas::Matrix{T}; power=2, rep=:polyhedron, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer, 
    T<:Real
}
    n = length(sample |> first)
    lib = default_library(n, T)

    return polyhedral_frechet_set(Opt, lib, sample, alphas; power=power, tol=tol)
end
function polyhedral_frechet_set(lib::Lib, sample, alphas; power=2, tol=1e-3) where {
    Lib<:Polyhedra.Library
}
    H = polyhedral_frechet_set(Clarabel.Optimizer, lib, sample, alphas; power=power, tol=tol) |> hrep
    redH = removehredundancy(H, Clarabel.Optimizer)
    return polyhedron(redH, lib)
end
polyhedral_frechet_set(sample, alphas; power=2, tol=1e-3) = polyhedral_frechet_set(CDDLib.Library(:exact), sample, alphas; power=power, tol=tol)


function polyhedral_ball(lib::Lib, center, alphas, radius::T) where {T<:Real, Lib<:Polyhedra.Library}
    b = radius .+ alphas * center

    return polyhedron(hrep(alphas, b), lib)
end


"""
    tropical_frechet_set(sample; power=2, rep=:polyhedron, tol=1e-3)

Find one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals scaled to α⋅x = 1.
`power` gives the exponent of the distance before taking the sum.
"""
function tropical_frechet_set(::Type{Opt}, lib::Lib, sample::Vector{Vector{T}}; power=2, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer, Lib<:Polyhedra.Library, T<:Real
}
    n = length(sample[1])
    alphas = tropical_ball_facets(n)

    P = polyhedral_frechet_set(Opt, lib, sample, alphas; power=power, tol=tol)
    return tropical_remove_redundant_halfspaces!(P)
end

function tropical_frechet_set(::Type{Opt}, sample::Vector{Vector{T}}; power=2, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer, T<:Real
}
    n = length(sample[1])
    lib = default_library(n, T)
    return tropical_frechet_set(Opt, lib, sample; power=power, tol=tol)
end
function tropical_frechet_set(lib::Lib, sample::Vector{Vector{T}}; power=2, tol=1e-3) where {
    Lib<:Polyhedra.Library, T<:Real
}
    H = tropical_frechet_set(Clarabel.Optimizer, lib, sample; power=power, tol=tol) |> hrep
    redH = removehredundancy(H, Clarabel.Optimizer)
    return polyhedron(redH, lib)
end
tropical_frechet_set(sample; power=2, tol=1e-3) = tropical_frechet_set(CDDLib.Library(:exact), sample; power=power, tol=tol)

function tropical_remove_redundant_halfspaces!(P::Polyhedron{T}) where T<:Real
    H = hrep(P)
    hs = [(h.a, h.β) for h in halfspaces(H)]
    hp = hyperplanes(H)

    bdict = Dict{Vector{T}, T}()
    
    for (a, β) in hs
        if a in keys(bdict)
            bdict[a] = min(bdict[a], β)
        else
            bdict[a] = β
        end
    end

    f = [HalfSpace(a, β) for (a, β) in bdict]
    Hnew = hrep(hp, f)

    return polyhedron(Hnew, library(P))
end

tropical_ball(lib::Lib, center::Vector{T}, radius::R) where {T<:Real, R<:Real, Lib<:Polyhedra.Library} =
    polyhedral_ball(lib, center, tropical_ball_facets(T, length(center)), radius)
tropical_ball(center::Vector{T}, radius::R) where {T<:Real, R<:Real} = tropical_ball(default_library(length(center), T), center, radius)

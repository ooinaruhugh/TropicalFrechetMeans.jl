using LinearAlgebra
using MathOptInterface, JuMP
using Polyhedra

function polyhedral_ball(center, alphas, radius::T) where {T<:Real}
    b = radius .+ alphas * center
  
    return alphas, b
end

polyhedral_ball(::Type{Polyhedron}, center, alphas, radius::T) where {T<:Real} = 
    hrep(polyhedral_ball(center, alphas, radius)...) |> polyhedron

polyhedral_ball(lib::Lib, center, alphas, radius::T) where {T<:Real, Lib<:Polyhedra.Library} =
    polyhedron(hrep(polyhedral_ball(center, alphas, radius)...), lib)

tropical_ball(lib::Lib, center::Vector{T}, radius::R) where {T<:Real, R<:Real, Lib<:Polyhedra.Library} =
    polyhedral_ball(lib, center, tropical_ball_facets(T, length(center)), radius)
tropical_ball(center::Vector{T}, radius::R) where {T<:Real, R<:Real} = tropical_ball(default_library(length(center), T), center, radius)


"""
Find the set of polyhedral Fréchet means of a given sample.
The rows of `alphas` are the facet normals scaled to α⋅x = 1.
`power` gives the exponent of the distance before taking the sum.
"""
function polyhedral_frechet_set(::Type{Opt}, lib::Lib, sample, alphas; power=2, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer,
    Lib<:Polyhedra.Library
}
    FM = polyhedral_frechet_mean(Opt, sample, alphas; power=power, tol=tol)
    @debug "Frechet mean found: $(FM)"

    b = foldl(sample; init=fill(Inf, size(alphas, 1))) do acc,pt
      r = polyhedral_distance(FM, pt, alphas)
      min.(acc, r .+ alphas * pt)
    end

    if tol > 0
      b = rationalize.(b; tol=tol)
    end
    
    return polyhedron(hrep(alphas, b), lib)
end
function polyhedral_frechet_set(::Type{Opt}, sample::Vector{Vector{T}}, alphas::Matrix{T}; power=2, rep=:polyhedron, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer, 
    T<:Real
}
    n = length(sample |> first)
    lib = default_library(n, T)

    return polyhedral_frechet_set(Opt, lib, sample, alphas; power=power, tol=tol)
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
  n = length(sample |> first)
  alphas = tropical_ball_facets(n)

  return polyhedral_frechet_set(Opt, lib, sample, alphas; power=power, tol=tol)
end

function tropical_frechet_set(::Type{Opt}, sample::Vector{Vector{T}}; power=2, tol=1e-3) where {
    Opt<:MathOptInterface.AbstractOptimizer, T<:Real
}
  n = length(sample |> first)
  lib = default_library(n, T)
  return tropical_frechet_set(Opt, lib, sample; power=power, tol=tol)
end

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

function breakpoints_of_tropical_line(p,q)
    r = q - p
    σ = sortperm(r)

    V = Vector{typeof(r)}(undef, length(r))

    for k in 1:length(r)
        V[k] = [q[σ[1:k]]..., r[σ[k]] .+ p[σ[k+1:end]]...][σ]
    end
    
    return V
end

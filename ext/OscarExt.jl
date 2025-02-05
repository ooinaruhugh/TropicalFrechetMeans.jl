module OscarExt

using TropicalFrechetMeans
using Polyhedra
import Oscar: QQ, TropicalSemiringElem

TropicalFrechetMeans.breakpoints_of_tropical_line(p::Vector{T}, q::Vector{T}) where {T<:TropicalSemiringElem} = breakpoints_of_tropical_line(QQ.(p), QQ.(q))

function TropicalFrechetMeans.wdp_from_polyhedron(H::Polyhedron)
    T = tropical_semiring()
    C = identity_matrix(T, n)
    for h in halfspaces(H)
        setindex!(C, h.β, indices(h.a)...)
    end

    return C
end


end # module OscarExt
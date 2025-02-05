module OscarExt

using TropicalFrechetMeans
import Oscar: QQ, TropicalSemiringElem

TropicalFrechetMeans.breakpoints_of_tropical_line(p::Vector{T}, q::Vector{T}) where {T<:TropicalSemiringElem} = breakpoints_of_tropical_line(QQ.(p), QQ.(q))


end # module OscarExt
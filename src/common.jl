
"""
    polyhedral_distance(vec1, vec2, alphas)

Calculate the polyhedral distance between two vectors.
The rows of `alphas` are the facet normals of the corresponding unit ball.
"""
function polyhedral_distance(vec1, vec2, alphas)
    differences = alphas * (vec1 - vec2)
    return maximum(differences)
end
# polyhedral_distance(x, y, P::Polyhedron) = polyhedral_distance(x, y, P.A)

"""
    tropical_distance(vec1, vec2)

Calculate the tropical distance between two vectors.
"""
function tropical_distance(vec1, vec2)
    return maximum(vec1 - vec2) - minimum(vec1 - vec2)
end

"""
    sum_of_poly_dist(ref, sample, alphas; power=1)

Calculate the sum of polyhedral distances between `ref` and the points in `sample`.
The rows of `alphas` are the facet normals of the corresponding unit ball.
`power` gives the exponent of the distance before taking the sum.
"""
function sum_of_poly_dist(ref, sample, alphas; power=1)
    return sum(polyhedral_distance(ref, s, alphas)^power for s in sample)
end

"""
    sum_of_trop_dist(ref, sample; power=1)

Calculate the sum of tropical distances between `ref` and the points in `sample`.
"""
function sum_of_trop_dist(ref, sample; power=1)
    return sum([tropical_distance(pt, ref)^power for pt in sample])
end

"""
    trop_facets(n::Int64)

Computes the facet normals in n dimensions for a tropical unit ball.
"""
function tropical_ball_facets(::Type{T}, n::Int64) where {T<:Real}
    result = zeros(T, n * (n - 1), n)
    k = 1
    for i = 1:n
        for j = 1:n
            if i != j
                result[k, i] = T(1)
                result[k, j] = T(-1)
                k += 1
            end
        end
    end
    return result
end
tropical_ball_facets(n::Int64) = tropical_ball_facets(Rational{Int64}, n)

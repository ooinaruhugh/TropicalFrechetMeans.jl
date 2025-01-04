module TropicalFrechetMeans

include("common.jl")
include("models.jl")
include("polyhedra.jl")

export polyhedral_distance
export sum_of_poly_dist
export tropical_distance
export sum_of_trop_dist

export trop_ball_facets

export polyhedral_frechet_model
export polyhedral_frechet_mean
export polyhedral_frechet_set

export tropical_frechet_model
export tropical_frechet_mean
export tropical_frechet_set

end # module TropicalFrechetMeans

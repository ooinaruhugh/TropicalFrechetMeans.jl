using MathOptInterface, JuMP, Clarabel

"""
    polyhedral_frechet_model(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}

Prepares a JuMP model for finding one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals of the unit ball for the polyhedral distance.
"""
function polyhedral_frechet_model(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}
    dim = length(sample[1])
    
    # Choose model depending on power
    if power == 1
        error("FW computation not yet implemeneted")
    end

    model = Model(Opt)
    
    # suppress printing
    set_silent(model)

    @variable(model, x[1:dim])
    @variable(model, t[1:length(sample)])

    @objective(model, Min, sum(t))

    for (p_idx, p) in enumerate(sample)
        expressions = alphas * (x - p)
        
        for expr in expressions
            @constraint(model, t[p_idx] >= (expr)^power)
        end
    end

    return model, x
end
polyhedral_frechet_model(sample, alphas; power=2) =polyhedral_frechet_model(Clarabel.Optimizer, sample, alphas; power=power)

"""
    polyhedral_frechet_mean(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}

Find one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals of the unit ball for the polyhedral distance.
`power` gives the exponent of the distance before taking the sum.
"""
function polyhedral_frechet_mean(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}
    model, x = polyhedral_frechet_model(Opt, sample, alphas, power=power)
    
    @debug "\nOptimising..."
    
    optimize!(model)
    minimiser = value.(x)
    
    return minimiser
end
polyhedral_frechet_mean(sample, alphas; power=2) = polyhedral_frechet_mean(Clarabel.Optimizer, sample, alphas; power=power)

"""
    tropical_frechet_model(::Type{Opt}, sample; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}

Prepares a JuMP model for finding one tropical Fréchet mean of a given sample.
`power` gives the exponent of the distance before taking the sum.
"""
function tropical_frechet_model(::Type{Opt}, sample; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}
    dim = length(sample[1])
    alphas = tropical_ball_facets(dim)
    return polyhedral_frechet_model(Opt, sample, alphas; power=power)
end

"""
    tropical_frechet_mean(::Type{Opt}, sample; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}

Find one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals scaled to α⋅x = 1.
`power` gives the exponent of the distance before taking the sum.
"""
function tropical_frechet_mean(::Type{Opt}, sample; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}
    model, x = tropical_frechet_model(Opt, sample; power=power)
    
    @debug "\nOptimising..."
    
    optimize!(model)
    minimiser = value.(x)
    
    return minimiser
end
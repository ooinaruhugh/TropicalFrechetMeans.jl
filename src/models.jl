using MathOptInterface, JuMP

"""
    polyhedral_frechet_model(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}

Prepares a JuMP model for finding one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals of the unit ball for the polyhedral distance.
"""
function polyhedral_frechet_model(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}
    dim = length(sample[1])
    
    # Choose model depending on power
    if power == 1
        error("Fermat-Weber computation not yet implemeneted")
    end

    model = Model(Opt)
    
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

"""
    polyhedral_frechet_mean(::Type{Opt}, sample, alphas; power=2) where {Opt<:MathOptInterface.AbstractOptimizer}

Find one polyhedral Fréchet mean of a given sample.
The rows of `alphas` are the facet normals of the unit ball for the polyhedral distance.
`power` gives the exponent of the distance before taking the sum.
"""
function polyhedral_frechet_mean(::Type{Opt}, sample, alphas; power=2, tol=1e-3) where {Opt<:MathOptInterface.AbstractOptimizer}
    if tol > 0 
        sample = map(sample) do pt; rationalize.(pt; tol=tol) end
    end
    model, x = polyhedral_frechet_model(Opt, sample, alphas, power=power)
    
    @debug "\nOptimising..."
    
    optimize!(model)
    minimiser = value.(x)
    
    if tol > 0 
        return rationalize.(minimiser; tol=tol)
    else 
        return minimiser
    end
end

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
function tropical_frechet_mean(::Type{Opt}, sample; power=2, tol=1e-3) where {Opt<:MathOptInterface.AbstractOptimizer}
    if tol > 0 
        sample = map(sample) do pt; rationalize.(pt; tol=tol) end
    end
    model, x = tropical_frechet_model(Opt, sample; power=power)
    
    @debug "\nOptimising..."
    
    optimize!(model)
    minimiser = value.(x)
    
    if tol > 0 
        return rationalize.(minimiser; tol=tol)
    else 
        return minimiser
    end
end

function conjugate_gradients(sample)
    alphas = tropical_ball_facets(length(sample[1])) |> transpose
    δ = x-> 2//(x+2)

    v = sum(sample) // length(sample)
    S = sum_of_trop_dist(v, sample; power=2)

    k=1
    while true
        g = alphas .+ δ(k)*v
        newS = [sum_of_trop_dist(g[:,k], sample; power=2) for k in axes(g, 2)]
        println(v, δ(k))
        println(g |> transpose)

        i = findfirst(newS .< S)
        if isnothing(i)
            break
        else
            S = newS[i]
            v = g[:,i]
            k += 1
        end
    end 

    return v
end

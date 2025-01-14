using LinearAlgebra
using Polyhedra
using TropicalFrechetMeans

using JSON3
using DataFrames

# Function to read the JSON file and convert to a list of matrices
function read_and_convert_json(file_path::String)
    # Read the JSON file
    json_data = JSON3.read(file_path)
    
    # Extract elements from the nested arrays
    elements = [x[1] for x in json_data]
    elements = rationalize.(10000 * elements, tol=1e-2)
    
    # Convert elements into matrices
    num_elements = length(elements)
    matrices = []
    
    for i in 1:64:num_elements
        # Get the next 64 elements
        matrix_elements = elements[i:min(i+63, num_elements)]
        
        # Convert to an 8x8 matrix if there are 64 elements, otherwise create a smaller matrix
        matrix_size = length(matrix_elements)
        sqrt_size = Int(sqrt(matrix_size))
        push!(matrices, reshape(matrix_elements, sqrt_size, sqrt_size))
    end
    
    return matrices
end

"""
Take a matrix of pairwise distances between taxa and returns the cophenetic vector.
"""
function cophenetic_from_distance(pairwise)
    n = size(pairwise, 1)
    coph = [pairwise[i, j] for i in 1:n-1 for j in i+1:n]
    return coph
end

# Read and convert the JSON file
file_path = "all_matrices.json"
matrices = read_and_convert_json(file_path)
taxa = ["Tg", "Et", "Cp", "Ta", "Bb", "Tt", "Pv", "Pf"]

coph_vecs = [cophenetic_from_distance(mat) for mat in matrices]

# @time @show tropical_frechet_mean(coph_vecs) 
@time phylo_frech = tropical_frechet_set(coph_vecs)
println(phylo_frech)
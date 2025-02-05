# Tropical Fréchet Means

TropicalFrechetMeans.jl is a Julia package providing algorithms to compute Tropical Fréchet means for a dataset.
This is the accompanying code to [link to arxiv].

## Usage

This package provides two main functions, numerically calculating a tropical Fréchet mean given a dataset, and calculating the
*Fréchet mean polytrope* given a dataset. For the numerical computation, we make use of the [JuMP ecosystem](https://jump.dev) and support 
quadratic solvers implementing the [MathOptInterface](https://jump.dev/MathOptInterface.jl/) (MOI).
The following shows a simple usage example and how to provide a MOI-compatible solver.

```julia
using Clarabel
using TropicalFrechetMeans

sample = [[-3,0,0], [0,-6,0], [0,0,-12]]

tropical_frechet_mean(Clarabel.Optimizer, sample)
```
The result of the above snippet should be `[0, 0, -1]`.

We recommend [Clarabel.jl](https://clarabel.org/) as a quadratic solver for use with this package. When loaded alongside this package, we
expose a default function that makes use of Clarabel. Thus, the computation of a tropical Fréchet mean is instead done by
```julia
tropical_frechet_mean(sample)
```

Calculating the entire Fréchet mean polytrope, one needs additionally an library for polyhedral computations implementing the interface of 
[Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl). One such example is the Julia version of [CDDLib](https://github.com/JuliaPolyhedra/CDDLib.jl).

Continuing above example, the corresponding Fréchet mean polytrope can be calculated in the following way.
```julia
using Clarabel, CDDLib
using TropicalFrechetMeans

P = tropical_frechet_set(sample)
```
Due to the nature of the computations, this will output a facet description of the Fréchet mean polytrope.
This can be converted to a vertex description, giving generators for all possible Fréchet means. 

## Contacts
The library is maintained by Kamillo Ferry (kafe (at) kafe (dot) dev).

## Acknowledgement
This implementation is based on code by Ariff Jazlan Johan.

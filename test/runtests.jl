using QuadraticOptimizer
using QuadraticOptimizer: M
using QuadraticOptimizer: center
using QuadraticOptimizer: hessian
using QuadraticOptimizer: distance
using StaticArrays
using LinearAlgebra
using Random
using Rotations
using ForwardDiff
using Test
using Aqua

Aqua.test_all(QuadraticOptimizer)

function F(p::StaticVector{D}) where D
    rng = Xoshiro(42)
    R = rand(rng, RotMatrix{D})
    A = R*Diagonal(3 .+ rand(rng, D))*R'
    b = randn(rng, D)
    C = randn(rng, D, D, D) / 20D
    return p'*A*p/2 + b'*p + sum(C[i,j,k]*p[i]*p[j]*p[k] for i in 1:D, j in 1:D, k in 1:D)
end
F(x::Real) = F(SVector(x))

function G(p::StaticVector{D}) where D
    rng = Xoshiro(42)
    A = Symmetric(rand(rng, -1:1, D, D) + 5I)
    b = rand(rng, -1:1, D)
    c = rand(rng, 0:5)
    q = Quadratic{D}(A,b,c)
    C = rand(rng, -1:1, D, D, D)
    return q(p) + sum(C[i,j,k]*p[i]*p[j]*p[k] for i in 1:D, j in 1:D, k in 1:D) / D
end
G(x::Real) = G(SVector(x))

function Q(::Val{D}) where D
    rng = Xoshiro(42)
    R = rand(rng, RotMatrix{D})
    A = Symmetric(R*Diagonal(1 .+ rand(rng, D))*R')
    b = randn(rng, D)
    return Quadratic{D}(A,b,0)
end

include("quadratic.jl")
include("qim.jl")
include("qfm.jl")
include("qim-qfm.jl")

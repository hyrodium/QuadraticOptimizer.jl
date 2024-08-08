module QuadraticOptimizer

using LinearAlgebra
using StaticArrays

export Quadratic
# export center
# export optimize!

include("quadratic.jl")
include("QIM_1-dim.jl")
include("QIM_n-dim.jl")

end

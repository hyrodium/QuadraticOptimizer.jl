module QuadraticOptimizer

using LinearAlgebra: issymmetric, tr, I, pinv
using StaticArrays: SVector, SHermitianCompact, SizedVector, SizedMatrix, StaticVector, StaticMatrix, SOneTo, arithmetic_closure

export Quadratic
export Linear
export quadratic_interpolation
export quadratic_fitting
export optimize_qim
export optimize_qfm
export optimize_qim!
export optimize_qfm!

include("quadratic.jl")
include("linear.jl")
include("QIM_1-dim.jl")
include("QFM_1-dim.jl")
include("QIM_n-dim.jl")
include("QFM_n-dim.jl")

end

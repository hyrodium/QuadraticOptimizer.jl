function _recursion_qfm!(xs::Vector{<:Real}, fs::Vector{<:Real}, X::Matrix, F::Vector)
    N = length(F)
    F[mod(length(xs), 1:N)] = fs[end]
    X[1,mod(length(xs), 1:N)] = xs[end]^2/2
    X[2,mod(length(xs), 1:N)] = xs[end]
    Y = (X*X')\(X*F)
    a = Y[1]
    b = Y[2]
    return -b/a
end

"""
    optimize_qfm!(f, xs::Vector{<:Real}, fs::Vector{<:Real}, n_iter::Integer) -> xs, fs

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `xs`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated. This vector will be updated in-place during the optimization process.
- `fs`: A vector of function values corresponding to the points in `xs`. This vector will be updated in-place during the optimization process.
- `n_iter`: The number of optimizing iterations. After execution, the length of `xs` will be `N + n_iter`, where `N = length(xs)` before execution.

!!! note
    In each step of the QFM, the last `N` points from `xs` and `fs` are used to interpolate with a quadratic function.
    The method iteratively refines the points and function values, extending `xs` and `fs` with `n_iter` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer

julia> f(x) = sin(x) + x^2/10  # Function to minimize
f (generic function with 1 method)

julia> xs_init = [1.2, 0.1, -2.2, -1.0]  # Initial points
4-element Vector{Float64}:
  1.2
  0.1
 -2.2
 -1.0

julia> xs = copy(xs_init);

julia> fs = f.(xs);

julia> optimize_qfm!(f, xs, fs, 5);  # Optimize 5 steps
```
"""
function optimize_qfm!(f, xs::Vector{T}, fs::Vector{T}, n_iter::Integer) where {T<:Real}
    U = arithmetic_closure(T)
    N = length(xs)
    length(xs) == N ≥ 3 || error("The length of initial values should be 3.")
    X = ones(U, 3, N)
    F = ones(U, N)
    X[1, 1:N-1] .= xs[1:N-1].^2/2
    X[2, 1:N-1] .= xs[1:N-1]
    F[1:N-1] .= fs[1:N-1]
    for _ in 1:n_iter
        x = _recursion_qfm!(xs, fs, X, F)
        push!(fs,f(x))
        push!(xs,x)
    end
    return xs, fs
end

"""
    optimize_qfm(f, xs_init::Vector{<:Real}, fs_init::Vector{<:Real}, n_iter::Integer) -> xs, fs

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `xs_init`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated.
- `fs_init`: A vector of function values corresponding to the points in `xs_init`.
- `n_iter`: The number of optimizing iterations. After execution, the length of `xs` will be `N + n_iter`, where `N = length(xs)` before execution.

!!! note
    In each step of the QFM, the last `N` points from `xs` and `fs` are used to interpolate with a quadratic function.
    The method iteratively refines the points and function values, extending `xs` and `fs` with `n_iter` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer

julia> f(x) = sin(x) + x^2/10  # Function to minimize
f (generic function with 1 method)

julia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)
3-element Vector{Float64}:
  1.2
  0.1
 -2.2

julia> xs, fs = optimize_qfm(f, xs_init, f.(xs_init), 5);  # Optimize 5 steps
```
"""
function optimize_qfm(f, xs_init::Vector{<:Real}, fs_init::Vector{<:Real}, n_iter::Integer)
    xs = copy(xs_init)
    fs = copy(fs_init)
    return optimize_qfm!(f, xs, fs, n_iter)
end

"""
    optimize_qfm(f, xs_init::Vector{<:Real}, n_iter::Integer) -> xs, fs

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `xs_init`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated.
- `n_iter`: The number of optimizing iterations. After execution, the length of `xs` will be `N + n_iter`, where `N = length(xs)` before execution.

!!! note
    In each step of the QFM, the last `N` points from `xs` and `fs` are used to interpolate with a quadratic function.
    The method iteratively refines the points and function values, extending `xs` and `fs` with `n_iter` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer

julia> f(x) = sin(x) + x^2/10  # Function to minimize
f (generic function with 1 method)

julia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)
3-element Vector{Float64}:
  1.2
  0.1
 -2.2

julia> xs = copy(xs_init)  # Keep initial points
3-element Vector{Float64}:
  1.2
  0.1
 -2.2

julia> xs, fs = optimize_qfm(f, xs, 5);  # Optimize 5 steps
```
"""
function optimize_qfm(f, xs_init::Vector{<:Real}, n_iter::Integer)
    xs = copy(xs_init)
    fs = f.(xs)
    return optimize_qfm!(f, xs, fs, n_iter)
end

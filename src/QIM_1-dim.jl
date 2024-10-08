function _recursion_qim!(xs::Vector{<:Real}, fs::Vector{<:Real}, X::StaticVector{3}, F::StaticVector{3})
    F[mod(length(xs), 1:3)] = fs[end]
    X[mod(length(xs), 1:3)] = xs[end]
    f₁ ,f₂, f₃ = F
    x₁, x₂, x₃ = X
    a₁ = f₁*(x₂-x₃)
    a₂ = f₂*(x₃-x₁)
    a₃ = f₃*(x₁-x₂)
    return (a₁*(x₂+x₃)+a₂*(x₃+x₁)+a₃*(x₁+x₂))/2(a₁+a₂+a₃)
end

"""
    optimize_qim!(f, xs::Vector{<:Real}, fs::Vector{<:Real}, n_iter::Integer) -> xs, fs

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `xs`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated. This vector will be updated in-place during the optimization process.
- `fs`: A vector of function values corresponding to the points in `xs`. This vector will be updated in-place during the optimization process.
- `n_iter`: The number of optimizing iterations. After execution, the length of `xs` will be `N + n_iter`, where `N = length(xs)` before execution.

!!! note
    In each step of the QIM, the last `3` points from `xs` and `fs` are used to interpolate with a quadratic function.
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

julia> xs = copy(xs_init);

julia> fs = f.(xs);

julia> optimize_qim!(f, xs, fs, 10);  # Optimize 10 steps
```
"""
function optimize_qim!(f, xs::Vector{T}, fs::Vector{T}, n_iter::Integer) where {T <: Real}
    length(xs) ≠ 3 && error("The length of initial values should be 3.")
    U = arithmetic_closure(T)
    X = SizedVector{3}(zeros(U, 3))
    F = SizedVector{3}(zeros(U, 3))
    X[1] = xs[1]
    X[2] = xs[2]
    F[1] = fs[1]
    F[2] = fs[2]
    for _ in 1:n_iter
        x = _recursion_qim!(xs, fs, X, F)
        push!(xs,x)
        push!(fs,f(x))
    end
    return xs, fs
end

"""
    optimize_qim(f, xs_init::Vector{<:Real}, fs_init::Vector{<:Real}, n_iter::Integer) -> xs, fs

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `xs_init`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated.
- `fs_init`: A vector of function values corresponding to the points in `xs_init`.
- `n_iter`: The number of optimizing iterations. After execution, the length of `xs` will be `N + n_iter`, where `N = length(xs)` before execution.

!!! note
    In each step of the QIM, the last `3` points from `xs` and `fs` are used to interpolate with a quadratic function.
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

julia> xs, fs = optimize_qim(f, xs_init, f.(xs_init), 10);  # Optimize 10 steps
```
"""
function optimize_qim(f, xs_init::Vector{<:Real}, fs_init::Vector{<:Real}, n_iter::Integer)
    xs = copy(xs_init)
    fs = copy(fs_init)
    return optimize_qim!(f, xs, fs, n_iter)
end

"""
    optimize_qim(f, xs_init::Vector{<:Real}, n_iter::Integer) -> xs, fs

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `xs_init`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated.
- `n_iter`: The number of optimizing iterations. After execution, the length of `xs` will be `N + n_iter`, where `N = length(xs)` before execution.

!!! note
    In each step of the QIM, the last `3` points from `xs` and `fs` are used to interpolate with a quadratic function.
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

julia> xs, fs = optimize_qim(f, xs, 10);  # Optimize 10 steps
```
"""
function optimize_qim(f, xs_init::Vector{<:Real}, n_iter::Integer)
    xs = copy(xs_init)
    fs = f.(xs)
    return optimize_qim!(f, xs, fs, n_iter)
end

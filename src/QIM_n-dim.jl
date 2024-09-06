"""
    quadratic_interpolation(ps::AbstractVector{<:StaticVector{D, <:Real}}, fs::AbstractVector{<:Real}) where D

Calculate quadratic interpolation based on input values.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> q = Quadratic{2}([2,1,3], [1,2], 2)
Quadratic{2, 3, Int64}([2, 1, 3], [1, 2], 2)

julia> ps = [@SVector rand(2) for _ in 1:6]
6-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]

julia> quadratic_interpolation(ps, q.(ps))
Quadratic{2, 3, Float64}([1.9999999999997884, 0.999999999999996, 2.999999999999987], [1.0000000000001212, 2.0000000000000053], 1.999999999999966)
```
"""
function quadratic_interpolation(ps::AbstractVector{<:StaticVector{D, T}}, fs::AbstractVector{T}) where {D, T<:Real}
    L = D*(D+1)÷2
    M = D+L+1
    N = length(ps)
    U = StaticArrays.arithmetic_closure(T)
    length(fs) == N == M || error("The length of initial values should be equal to $(M).")
    X = SizedMatrix{M,M}(ones(U, M, M))
    F = SizedVector{M}(zeros(U, M))
    _initialize_XF!(X, F, ps, fs)
    _update_XF_at_j!(X, F, ps, fs, mod(length(ps), 1:M))
    return _quadratic(X, F, Val(D))
end

Base.@deprecate interpolation quadratic_interpolation false

"""
    optimize_qim!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `ps`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated. This vector will be updated in-place during the optimization process.
- `fs`: A vector of function values corresponding to the points in `ps`. This vector will be updated in-place during the optimization process.
- `n`: The number of optimizing iterations. After execution, the length of `ps` will be `m + n`, where `m = length(ps)` before execution.

!!! note
    In each step of the QIM, the last `M` (`==((D+2)*(D+1)/2)`) points from `ps` and `fs` are used to interpolate with a quadratic function.
    The method iteratively refines the points and function values, extending `ps` and `fs` with `n` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
f (generic function with 1 method)

julia> ps_init = [@SVector rand(2) for _ in 1:6]
6-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]

julia> ps = copy(ps_init);

julia> fs = f.(ps);

julia> optimize_qim!(f, ps, fs, 20);
```
"""
function optimize_qim!(f, ps::Vector{<:SVector{D, T}}, fs::Vector{T}, n::Integer) where {D, T<:Real}
    L = D*(D+1)÷2
    M = D+L+1
    N = length(ps)
    U = StaticArrays.arithmetic_closure(T)
    length(fs) == N == M || error("The length of initial values should be equal to $(M).")
    X = SizedMatrix{M,M}(ones(U, M, M))
    F = SizedVector{M}(zeros(U, M))
    _initialize_XF!(X, F, ps, fs)
    for _ in 1:n
        _update_XF_at_j!(X, F, ps, fs, mod(length(ps), 1:M))
        q = _quadratic(X, F, Val(D))
        p = center(q)
        push!(fs,f(p))
        push!(ps,p)
    end
    return ps, fs
end

"""
    optimize_qim(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `ps`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated.
- `fs`: A vector of function values corresponding to the points in `ps`.
- `n`: The number of optimizing iterations. After execution, the length of `ps` will be `m + n`, where `m = length(ps)` before execution.

!!! note
    In each step of the QIM, the last `M` (`==((D+2)*(D+1)/2)`) points from `ps` and `fs` are used to interpolate with a quadratic function.
    The method iteratively refines the points and function values, extending `ps` and `fs` with `n` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
f (generic function with 1 method)

julia> ps_init = [@SVector rand(2) for _ in 1:6]
6-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]

julia> ps, fs = optimize_qim(f, ps_init, f.(ps_init), 20);
```
"""
function optimize_qim(f, ps_init::Vector{<:SVector{D, <:Real}}, fs_init::Vector{<:Real}, n::Integer) where D
    ps = copy(ps_init)
    fs = copy(fs_init)
    return optimize_qim!(f, ps, fs, n)
end

"""
    optimize_qim(f, ps::Vector{<:SVector{D, <:Real}}, n::Integer)

Optimize a function `f` using the Quadratic Interpolation Method (QIM).

# Arguments
- `f`: The objective function to be optimized.
- `ps`: A vector of points in ``\\mathbb{R}^D`` where `f` has been evaluated.
- `n`: The number of optimizing iterations. After execution, the length of `ps` will be `m + n`, where `m = length(ps)` before execution.

!!! note
    In each step of the QIM, the last `M` (`==((D+2)*(D+1)/2)`) points from `ps` and `fs` are used to interpolate with a quadratic function.
    The method iteratively refines the points and function values, extending `ps` and `fs` with `n` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
f (generic function with 1 method)

julia> ps_init = [@SVector rand(2) for _ in 1:6]
6-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]

julia> ps, fs = optimize_qim(f, ps_init, 10);
```
"""
function optimize_qim(f, ps_init::Vector{<:SVector{D, <:Real}}, n::Integer) where D
    ps = copy(ps_init)
    fs = f.(ps)
    return optimize_qim!(f, ps, fs, n)
end

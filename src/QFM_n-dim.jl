"""
    quadratic_fitting(ps::AbstractVector{<:StaticVector{D, <:Real}}, fs::AbstractVector{<:Real}) where D

Calculate quadratic fitting based on input values.

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

julia> quadratic_fitting(vcat(ps, ps), vcat(q.(ps).-1, q.(ps).+1))
Quadratic{2, 3, Float64}([2.000000000387182, 0.9999999999890256, 2.999999999972893], [0.9999999997884512, 2.0000000000207985], 2.000000000052808)
```
"""
function quadratic_fitting(ps::AbstractVector{<:StaticVector{D, T}}, fs::AbstractVector{T}) where {D, T<:Real}
    L = D*(D+1)÷2
    M = D+L+1
    N = length(ps)
    U = StaticArrays.arithmetic_closure(T)
    length(fs) == N ≥ M || error("The length of initial values should be larger than $(M).")
    X = ones(U, M, N)
    F = zeros(U, N)
    _initialize_XF!(X, F, ps, fs)
    _update_XF_at_j!(X, F, ps, fs, mod(length(ps), 1:N))
    return _quadratic(X, F, Val(D))
end

Base.@deprecate fitting quadratic_fitting false

"""
    optimize_qfm!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)

Optimize a function `f` using the Quadratic Fitting Method (QFM).

# Arguments
- `f`: The objective function to be optimized.
- `ps`: A vector of points in ``\\mathbb{R}^d`` where `f` has been evaluated. This vector will be updated in-place during the optimization process.
- `fs`: A vector of function values corresponding to the points in `ps`. This vector will be updated in-place during the optimization process.
- `n`: The number of optimizing iterations. After execution, the length of `ps` will be `m + n`, where `m = length(ps)` before execution.

!!! note
    In each step of the QFM, the last `m` points from `ps` and `fs` are used to fit a quadratic function.
    The method iteratively refines the points and function values, extending `ps` and `fs` with `n` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
f (generic function with 1 method)

julia> ps_init = [@SVector rand(2) for _ in 1:10]
10-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]
 [0.34264792290134793, 0.2678704383989201]
 [0.515871408349502, 0.09002301604339691]
 [0.27265744579429385, 0.191562202596938]
 [0.4235912564725989, 0.4847023673932017]

julia> ps = copy(ps_init);

julia> fs = f.(ps);

julia> optimize_qfm!(f, ps, fs, 20);
```
"""
function optimize_qfm!(f, ps::Vector{<:SVector{D, T}}, fs::Vector{T}, n::Integer) where {D, T<:Real}
    L = D*(D+1)÷2
    M = D+L+1
    N = length(ps)
    U = StaticArrays.arithmetic_closure(T)
    length(fs) == N ≥ M || error("The length of initial values should be larger than $(M).")
    X = ones(U, M, N)
    F = zeros(U, N)
    _initialize_XF!(X, F, ps, fs)
    for _ in 1:n
        _update_XF_at_j!(X, F, ps, fs, mod(length(ps), 1:N))
        q = _quadratic(X, F, Val(D))
        p = center(q)
        push!(fs,f(p))
        push!(ps,p)
    end
    return ps, fs
end

"""
    optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)

Optimize a function `f` using the Quadratic Fitting Method (QFM).

# Arguments
- `f`: The objective function to be optimized.
- `ps`: A vector of points in ``\\mathbb{R}^d`` where `f` has been evaluated.
- `fs`: A vector of function values corresponding to the points in `ps`.
- `n`: The number of optimizing iterations. After execution, the length of `ps` will be `m + n`, where `m = length(ps)` before execution.

!!! note
    In each step of the QFM, the last `m` points from `ps` and `fs` are used to fit a quadratic function.
    The method iteratively refines the points and function values, extending `ps` and `fs` with `n` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
f (generic function with 1 method)

julia> ps_init = [@SVector rand(2) for _ in 1:10]
10-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]
 [0.34264792290134793, 0.2678704383989201]
 [0.515871408349502, 0.09002301604339691]
 [0.27265744579429385, 0.191562202596938]
 [0.4235912564725989, 0.4847023673932017]

julia> ps, fs = optimize_qfm(f, ps_init, f.(ps_init), 20);
```
"""
function optimize_qfm(f, ps_init::Vector{<:SVector{D, <:Real}}, fs_init::Vector{<:Real}, n::Integer) where D
    ps = copy(ps_init)
    fs = copy(fs_init)
    return optimize_qfm!(f, ps, fs, n)
end

"""
    optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, n::Integer)

Optimize a function `f` using the Quadratic Fitting Method (QFM).

# Arguments
- `f`: The objective function to be optimized.
- `ps`: A vector of points in ``\\mathbb{R}^d`` where `f` has been evaluated.
- `n`: The number of optimizing iterations. After execution, the length of `ps` will be `m + n`, where `m = length(ps)` before execution.

!!! note
    In each step of the QFM, the last `m` points from `ps` and `fs` are used to fit a quadratic function.
    The method iteratively refines the points and function values, extending `ps` and `fs` with `n` additional points resulting from the optimization process.

# Examples
```jldoctest
julia> using QuadraticOptimizer, StaticArrays, Random

julia> Random.seed!(42);

julia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
f (generic function with 1 method)

julia> ps_init = [@SVector rand(2) for _ in 1:10]
10-element Vector{SVector{2, Float64}}:
 [0.6293451231426089, 0.4503389405961936]
 [0.47740714343281776, 0.7031298490032014]
 [0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278]
 [0.4570310908017041, 0.2993652953937611]
 [0.6611433726193705, 0.6394313620423493]
 [0.34264792290134793, 0.2678704383989201]
 [0.515871408349502, 0.09002301604339691]
 [0.27265744579429385, 0.191562202596938]
 [0.4235912564725989, 0.4847023673932017]

julia> ps, fs = optimize_qfm(f, ps_init, 20);
```
"""
function optimize_qfm(f, ps_init::Vector{<:SVector{D, <:Real}}, n::Integer) where D
    ps = copy(ps_init)
    fs = f.(ps)
    return optimize_qfm!(f, ps, fs, n)
end

"""
    Quadratic{D, T<:Real, L}
    Quadratic(a,b,c)

A struct that represents quadratic polynomial on ``\\mathbb{R}^D``.

```math
\\begin{aligned}
q(p) &= \\frac{1}{2}p^{\\intercal} A p + b^{\\intercal} p + c \\\\
A &= \\begin{pmatrix}
a_1     & a_2      & \\cdots & a_D \\\\
a_2     & a_{D+1}  & \\cdots & a_{2D-1} \\\\
\\vdots & \\vdots  & \\ddots & \\vdots \\\\
a_D     & a_{2D-1} & \\cdots & a_{L}
\\end{pmatrix}
\\end{aligned}
```

# Examples
```jldoctest
julia> q = Quadratic{2}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)
Quadratic{2, 3, Float64}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)

julia> q([1,2])
7.7
```
"""
struct Quadratic{D, T<:Real, L}
    a::SVector{L, T}
    b::SVector{D, T}
    c::T
    function Quadratic{D,T,L}(a::SVector{L, T}, b::SVector{D, T}, c::T) where {D, T<:Real, L}
        L ≠ D*(D+1)÷2 && throw(ArgumentError("The sizes of input vectors are invalid"))
        return new{D,T,L}(a, b, c)
    end
end

# (a,b,c) constructor
function Quadratic(a::StaticVector{L,Ta}, b::StaticVector{D,Tb}, c::Tc) where {L, D, Ta<:Real, Tb<:Real, Tc<:Real}
    T = promote_type(Ta, Tb, Tc)
    return Quadratic{D,T,L}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end
function Quadratic{D}(a::AbstractVector{Ta}, b::AbstractVector{Tb}, c::Tc) where {D, Ta<:Real, Tb<:Real, Tc<:Real}
    T = promote_type(Ta, Tb, Tc)
    L = D*(D+1)÷2
    return Quadratic{D,T,L}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end
function Quadratic{D,T}(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}, c::Real) where {D, T<:Real}
    L = D*(D+1)÷2
    return Quadratic{D,T,L}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end
function Quadratic{D,T,L}(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}, c::Real) where {D, T, L}
    return Quadratic{D,T,L}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end

# (A,b,c) constructor
function (::Type{Q})(A::AbstractMatrix{Ta}, b::AbstractVector{Tb}, c::Tc) where {Q<:Quadratic{D}, Ta<:Real, Tb<:Real, Tc<:Real} where D
    issymmetric(A) || throw(ArgumentError("Input matrix must be symmetric."))
    a = SVector(SHermitianCompact{D}(A).lowertriangle)
    return Q(a, b, c)
end
function Quadratic(A::StaticMatrix{D,D,Ta}, b::StaticVector{D,Tb}, c::Tc) where {Ta<:Real, Tb<:Real, Tc<:Real} where D
    issymmetric(A) || throw(ArgumentError("Input matrix must be symmetric."))
    a = SVector(SHermitianCompact{D}(A).lowertriangle)
    return Quadratic(a, b, c)
end
# Not sure the following methods causes ambiguities..
# function Quadratic(A::StaticMatrix{D,D,Ta}, b::AbstractVector{Tb}, c::Tc) where {Ta<:Real, Tb<:Real, Tc<:Real} where D
#     issymmetric(A) || throw(ArgumentError("Input matrix must be symmetric."))
#     a = SVector(SHermitianCompact{D}(A).lowertriangle)
#     return Quadratic(a, b, c)
# end
# function Quadratic(A::AbstractMatrix{Ta}, b::StaticVector{D,Tb}, c::Tc) where {Ta<:Real, Tb<:Real, Tc<:Real} where D
#     issymmetric(A) || throw(ArgumentError("Input matrix must be symmetric."))
#     a = SVector(SHermitianCompact{D}(A).lowertriangle)
#     return Quadratic(a, b, c)
# end

# (c) contructor
function Quadratic{D,T,L}(c) where {D, T<:Real, L}
    return Quadratic{D,T,L}(zero(SVector{L,T}), zero(SVector{D,T}), T(c))
end
function Quadratic{D,T}(c) where {D, T<:Real}
    L = D*(D+1)÷2
    return Quadratic{D,T,L}(zero(SVector{L,T}), zero(SVector{D,T}), T(c))
end
function Quadratic{D}(c::T) where {D, T<:Real}
    L = D*(D+1)÷2
    return Quadratic{D,T,L}(zero(SVector{L,T}), zero(SVector{D,T}), c)
end

Base.:(==)(q1::Quadratic, q2::Quadratic) = (q1.a == q2.a) & (q1.b == q2.b) & (q1.c == q2.c)

function Base.convert(::Type{Quadratic{D,T,L}}, c::Real) where {D, T<:Real, L}
    return Quadratic{D,T,L}(c)
end

function Base.promote_rule(::Type{Quadratic{D,T,L}}, ::Type{S}) where {D,T,L,S}
    return Quadratic{D,promote_type(T,S),L}
end

function Base.iszero(q::Quadratic)
    return iszero(q.a) & iszero(q.b) & iszero(q.c)
end

function Base.isfinite(q::Quadratic)
    return all(isfinite.(q.a)) & all(isfinite.(q.b)) & isfinite(q.c)
end

Base.:+(q::Quadratic) = q
Base.:-(q::Quadratic) = Quadratic(-q.a, -q.b, -q.c)
Base.:+(q1::Quadratic, q2::Quadratic) = Quadratic(q1.a+q2.a, q1.b+q2.b, q1.c+q2.c)
Base.:-(q1::Quadratic, q2::Quadratic) = Quadratic(q1.a-q2.a, q1.b-q2.b, q1.c-q2.c)
Base.:*(k::Real, q::Quadratic) = Quadratic(k*q.a, k*q.b, k*q.c)
Base.:*(q::Quadratic, k::Real) = k*q
Base.:\(k::Real, q::Quadratic) = Quadratic(k\q.a, k\q.b, k\q.c)
Base.:/(q::Quadratic, k::Real) = Quadratic(q.a/k, q.b/k, q.c/k)
Base.zero(::Type{Quadratic{D,T,L}}) where {D,T,L} = Quadratic{D,T,L}(zero(T))
Base.zero(::Type{Quadratic{D,T}}) where {D,T} = Quadratic{D,T}(zero(T))
Base.zero(::Type{Quadratic{D}}) where D = Quadratic{D,Float64}(zero(Float64))
Base.zero(q::Quadratic) = zero(typeof(q))

norm_quadratic(q::Quadratic) = sqrt(tr(hessian(q)*hessian(q))+q(center(q))^2)

function Base.rtoldefault(::Union{Q1,Type{Q1}}, ::Union{Q2,Type{Q2}}, atol::Real) where {Q1<:Quadratic{D,T1,L},Q2<:Quadratic{D,T2,L}} where {D,L,T1,T2}
    rtol = max(Base.rtoldefault(T1), Base.rtoldefault(T2))
    return atol > 0 ? zero(rtol) : rtol
end

function Base.isapprox(q1::Quadratic{D}, q2::Quadratic{D};
    atol::Real=0, rtol::Real=Base.rtoldefault(q1,q2,atol),
    nans::Bool=false, norm::Function=norm_quadratic) where D
    q1 == q2 ||
    (isfinite(q1) && isfinite(q2) && norm(q1-q2) <= max(atol, rtol*max(norm(q1), norm(q2)))) ||
    (nans && isnan(q1) && isnan(q2))
end

function center(q::Quadratic)
    return -hessian(q)\q.b
end

"""
    hessian(q::Quadratic)

Calculate the (constant) hessian of the quadratic polynomial.

# Examples
```jldoctest
julia> using QuadraticOptimizer: hessian

julia> using ForwardDiff

julia> q = Quadratic{2}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)
Quadratic{2, 3, Float64}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)

julia> hessian(q)
2×2 StaticArrays.SHermitianCompact{2, Float64, 3} with indices SOneTo(2)×SOneTo(2):
 2.0  1.0
 1.0  3.0

julia> ForwardDiff.hessian(q, [1,2])
2×2 Matrix{Float64}:
 2.0  1.0
 1.0  3.0
```
"""
function hessian(q::Quadratic)
    return SHermitianCompact(q.a)
end

function (q::Quadratic)(p)
    b = q.b
    c = q.c
    A = hessian(q)
    return p'*A*p/2 + b'*p + c
end

function _update_XF_at_j!(X::AbstractMatrix, F::AbstractVector, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, j::Integer) where {D}
    L = D*(D+1)÷2
    p = ps[end]
    F[j] = fs[end]
    i = 1
    for i1 in 1:D, i2 in i1:D
        if i1 == i2
            X[i,j] = p[i1]*p[i2]/2
        else
            X[i,j] = p[i1]*p[i2]
        end
        i = i + 1
    end
    X[L+1:L+D, j] .= p
    return X, F
end

function _initialize_XF!(X::AbstractMatrix, F::AbstractVector, ps::Vector{<:SVector{D, <:Real}}, fs) where {D}
    L = D*(D+1)÷2
    M, N = size(X)
    for j in 1:N-1
        p = ps[j]
        F[j] = fs[j]
        i = 1
        for i1 in 1:D, i2 in i1:D
            if i1 == i2
                X[i,j] = p[i1]*p[i2]/2
            else
                X[i,j] = p[i1]*p[i2]
            end
            i = i + 1
        end
        X[L+1:L+D, j] .= p
    end
    return X, F
end

"""
    M(D::Integer) = (D+2)*(D+1)÷2

Calculate the number of required points for QIM/QFM that is equal to the number of terms in `D`-dimensional quadratic polynomial.
```math
\\begin{aligned}
M(1) &= 3 \\quad (\\text{Number of terms in} \\ \\frac{a_1}{2} x^2 + b_1 x + c) \\\\
M(2) &= 6 \\quad (\\text{Number of terms in} \\ \\frac{a_1}{2} x^2 + a_2 xy + \\frac{a_3}{2} y^2 + b_1 x + b_2 y + c) \\\\
M(0) &= 1 \\quad (\\text{Number of terms in} \\ c)
\\end{aligned}
```
# Examples
```jldoctest
julia> using QuadraticOptimizer: M

julia> M(1)  # 3 points to interpolate a parabola
3

julia> M(2)  # 6 points to interpolate a 2-dimensional quadratic
6

julia> M(0)  # 1 point for constant, terms
1
```
"""
M(D::Integer) = (D+2)*(D+1)÷2

function _quadratic(X::StaticMatrix{N,N}, F::StaticVector{N}, ::Val{D}) where {N, D}
    N ≠ M(D) && throw(ArgumentError("The input value N=$(N) must be equal to M=(D+2)(D+1)/2=$(M(D))."))
    L = D*(D+1)÷2
    Y = X' \ F
    a = Y[SOneTo(L)]
    b = SVector{D}(Y[L+1:L+D])
    c = Y[end]
    return Quadratic(a,b,c)
end

function _quadratic(X::AbstractMatrix, F::AbstractVector, ::Val{D}) where {D}
    L = D*(D+1)÷2
    Y = (X*X')\(X*F)
    a = Y[SOneTo(L)]
    b = SVector{D}(Y[L+1:L+D])
    c = Y[end]
    return Quadratic(a,b,c)
end

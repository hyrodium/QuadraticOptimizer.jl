"""
    Linear{D, T<:Real}
    Linear(a, b)

A struct that represents linear polynomial on ``\\mathbb{R}^D``.

```math
l(p) = a^{\\intercal} p + b
```

# Examples
```jldoctest
julia> l = Linear{2}([-1.2, -2.3], 4.5)
Linear{2, Float64}([-1.2, -2.3], 4.5)

julia> l([1, 2])
-1.3
```
"""
struct Linear{D, T<:Real}
    a::SVector{D, T}
    b::T
    function Linear{D,T}(a::SVector{D, T}, b::T) where {D, T<:Real}
        return new{D,T}(a, b)
    end
end

# (a, b) constructor
function Linear(a::StaticVector{D,Ta}, b::Tb) where {D, Ta<:Real, Tb<:Real}
    T = promote_type(Ta, Tb)
    return Linear{D,T}(SVector{D,T}(a), T(b))
end
function Linear{D}(a::AbstractVector{Ta}, b::Tb) where {D, Ta<:Real, Tb<:Real}
    T = promote_type(Ta, Tb)
    return Linear{D,T}(SVector{D,T}(a), T(b))
end
function Linear{D,T}(a::AbstractVector{<:Real}, b::Real) where {D, T<:Real}
    return Linear{D,T}(SVector{D,T}(a), T(b))
end

# (b) constructor (constant only, a = 0)
function Linear{D,T}(b) where {D, T<:Real}
    return Linear{D,T}(zero(SVector{D,T}), T(b))
end
function Linear{D}(b::T) where {D, T<:Real}
    return Linear{D,T}(zero(SVector{D,T}), b)
end

Base.:(==)(l1::Linear, l2::Linear) = (l1.a == l2.a) & (l1.b == l2.b)

function Base.promote_rule(::Type{Linear{D,T}}, ::Type{S}) where {D,T,S}
    return Linear{D,promote_type(T,S)}
end

function Base.iszero(l::Linear)
    return iszero(l.a) & iszero(l.b)
end

function Base.isfinite(l::Linear)
    return all(isfinite.(l.a)) & isfinite(l.b)
end

function Base.isnan(l::Linear)
    return any(isnan.(l.a)) | isnan(l.b)
end

Base.:+(l::Linear) = l
Base.:-(l::Linear) = Linear(-l.a, -l.b)
Base.:+(l1::Linear, l2::Linear) = Linear(l1.a+l2.a, l1.b+l2.b)
Base.:-(l1::Linear, l2::Linear) = Linear(l1.a-l2.a, l1.b-l2.b)
Base.:+(l::Linear, c::Real) = Linear(l.a, l.b + c)
Base.:-(l::Linear, c::Real) = Linear(l.a, l.b - c)
Base.:+(c::Real, l::Linear) = Linear(l.a, l.b + c)
Base.:-(c::Real, l::Linear) = Linear(-l.a, c - l.b)
Base.:*(k::Real, l::Linear) = Linear(k*l.a, k*l.b)
Base.:*(l::Linear, k::Real) = k*l
Base.:\(k::Real, l::Linear) = Linear(k\l.a, k\l.b)
Base.:/(l::Linear, k::Real) = Linear(l.a/k, l.b/k)
Base.zero(::Type{Linear{D,T}}) where {D,T} = Linear{D,T}(zero(T))
Base.zero(::Type{Linear{D}}) where D = Linear{D,Float64}(zero(Float64))
Base.zero(l::Linear) = zero(typeof(l))

function (l::Linear)(p)
    return l.a'*p + l.b
end

# TODO: distance(l1::Linear, l2::Linear)
# TODO: Base.rtoldefault for Linear
# TODO: Base.isapprox for Linear

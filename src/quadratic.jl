struct Quadratic{D, L, T<:Real}
    a::SVector{L, T}
    b::SVector{D, T}
    c::T
end

function Quadratic(a::StaticVector{L,Ta}, b::StaticVector{D,Tb}, c::Tc) where {L, D, Ta<:Real, Tb<:Real, Tc <: Real}
    T = promote_type(Ta, Tb, Tc)
    Quadratic{D,L,T}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end

function Quadratic{D}(a::AbstractVector{Ta}, b::AbstractVector{Tb}, c::Tc) where {D, Ta<:Real, Tb<:Real, Tc <: Real}
    T = promote_type(Ta, Tb, Tc)
    L = D*(D+1)÷2
    Quadratic{D,L,T}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end

function Quadratic{D,L}(a::AbstractVector{Ta}, b::AbstractVector{Tb}, c::Tc) where {D, L, Ta<:Real, Tb<:Real, Tc <: Real}
    T = promote_type(Ta, Tb, Tc)
    Quadratic{D,L,T}(SVector{L,T}(a), SVector{D,T}(b), T(c))
end

function Base.:≈(q1::Quadratic, q2::Quadratic)
    (q1.a ≈ q2.a) & (q1.b ≈ q2.b) & (q1.c ≈ q2.c)
end

function Base.iszero(q::Quadratic)
    return iszero(q.a) & iszero(q.b) & iszero(q.c)
end

function Base.isfinite(q::Quadratic)
    return all(isfinite.(q.a)) & all(isfinite.(q.b)) & isfinite(q.c)
end

function Base.:+(q1::Quadratic, q2::Quadratic)
    return Quadratic(q1.a+q2.a, q1.b+q2.b, q1.c+q2.c)
end

function Base.:-(q1::Quadratic, q2::Quadratic)
    return Quadratic(q1.a-q2.a, q1.b-q2.b, q1.c-q2.c)
end

function Base.:+(q::Quadratic)
    return q
end

function Base.:-(q::Quadratic)
    return Quadratic(-q.a, -q.b, -q.c)
end

function Base.:*(k::Real, q::Quadratic)
    return Quadratic(k*q.a, k*q.b, k*q.c)
end

function Base.:*(q::Quadratic, k::Real)
    return k*q
end

function Base.:\(k::Real, q::Quadratic)
    return Quadratic(k\q.a, k\q.b, k\q.c)
end

function Base.:/(q::Quadratic, k::Real)
    return Quadratic(q.a/k, q.b/k, q.c/k)
end

# TODO
# function isapprox(q1::Quadratic{L,T1}, q2::Quadratic{L,T2};
#     atol::Real=0, rtol::Real=rtoldefault(q1,q2,atol),
#     nans::Bool=false, norm::Function=abs)
#     q1 == q2 ||
#     (isfinite(q1) && isfinite(q2) && norm(q1-q2) <= max(atol, rtol*max(norm(q1), norm(q2)))) ||
#     (nans && isnan(q1) && isnan(q2))
# end

function center(q::Quadratic{D}) where {D}
    a = q.a
    b = q.b

    A = MMatrix{D,D}(zeros(D,D))
    j = 1
    for i1 in 1:D
        for i2 in i1:D
            if i1 == i2
                A[i1,i2] = 2a[j]
            else
                A[i1,i2] = A[i2,i1] = a[j]
            end
            j = j + 1
        end
    end
    return -A\b
end

function (q::Quadratic{D,L,T})(p) where {D,L,T}
    a = q.a
    b = q.b
    c = q.c
    y = zero(T)
    j = 1
    for i1 in 1:D, i2 in i1:D
        y += a[j]*p[i1]*p[i2]
        j = j + 1
    end
    for i in 1:D
        y += b[i]*p[i]
    end
    y += c
    return y
end

function _update_XF_at_j!(X::AbstractMatrix, F::AbstractVector, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, j::Integer) where {D}
    L = D*(D+1)÷2
    p = ps[end]
    F[j] = fs[end]
    i = 1
    for i1 in 1:D, i2 in i1:D
        X[i,j] = p[i1]*p[i2]
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
            X[i,j] = p[i1]*p[i2]
            i = i + 1
        end
        X[L+1:L+D, j] .= p
    end
    return X, F
end

function _quadratic(X::StaticMatrix{N,N}, F::StaticVector{N}, ::Val{D}) where {N, D}
    L = D*(D+1)÷2
    Y = pinv(X') * F
    a = Y[SOneTo(L)]
    b = SVector{D}(Y[L+1:L+D])
    c = Y[end]
    return Quadratic(a,b,c)
end

function _quadratic(X::AbstractMatrix, F::AbstractVector, ::Val{D}) where {D}
    L = D*(D+1)÷2
    Y = pinv(X*X')*(X*F)
    a = Y[SOneTo(L)]
    b = SVector{D}(Y[L+1:L+D])
    c = Y[end]
    return Quadratic(a,b,c)
end

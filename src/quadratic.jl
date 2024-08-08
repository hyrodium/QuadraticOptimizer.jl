struct Quadratic{D, M, T<:Real}
    a::SVector{M, T}
    b::SVector{D, T}
    c::T
end

function Quadratic(a::StaticVector{M,Ta},b::StaticVector{D, Tb},c::Tc) where {M, D, Ta<:Real, Tb<:Real, Tc <: Real}
    T = promote_type(Ta, Tb, Tc)
    Quadratic{D,M,T}(SVector{M,T}(a), SVector{D,T}(b), T(c))
end

function center(q::Quadratic{D,M}) where {D,M}
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

function (q::Quadratic{D,M,T})(p) where {D,M,T}
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

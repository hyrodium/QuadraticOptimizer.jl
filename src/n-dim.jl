function _recursion!(f, ps::Vector{<:SVector{D, <:Real}}, F::MVector{N}, X::MMatrix{N,N}) where {D, N}
    M = D*(D+1)÷2
    p = ps[end]
    i = mod(length(ps), 1:N)
    F[i] = f(p...)
    j = 1
    for i1 in 1:D, i2 in i1:D
        X[i,j] = p[i1]*p[i2]
        j = j + 1
    end
    X[i,M+1:M+D] .= p
    Y = pinv(X) * F
    a = Y[1:M]
    b = SVector{D}(Y[M+1:M+D])
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

function optimize!(f, ps::Vector{<:SVector{D, <:Real}}, n::Integer) where D
    M = D*(D+1)÷2
    N = D+M+1
    length(ps) ≠ N && error("The length of initial values should be $N.")
    F = @MVector zeros(N)
    X = @MMatrix ones(N,N)
    for i in 1:N-1
        p = ps[i]
        F[i] = f(p...)
        j = 1
        for i1 in 1:D, i2 in i1:D
            X[i,j] = p[i1]*p[i2]
            j = j + 1
        end
        X[i,M+1:M+D] .= p
        X[i,N] = 1
    end
    for _ in 1:n
        p = _recursion!(f, ps, F, X)
        push!(ps,p)
    end
    return ps
end

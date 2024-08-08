function _recursion!(f, ps::Vector{<:SVector{D, <:Real}}, F::Vector, X::Matrix, L::Integer) where {D, N}
    M = D*(D+1)÷2
    p = ps[end]
    i = mod(length(ps), 1:L)
    F[i] = f(p)
    j = 1
    for i1 in 1:D, i2 in i1:D
        X[i,j] = p[i1]*p[i2]
        j = j + 1
    end
    X[i,M+1:M+D] .= p
    Y = pinv(X'*X)*(X'*F)
    a = Y[SOneTo(M)]
    b = SVector{D}(Y[M+1:M+D])
    c = Y[end]
    q = Quadratic(a,b,c)
    return center(q)
end

function optimize_qfm!(f, ps::Vector{<:SVector{D, <:Real}}, n::Integer) where D
    L = length(ps)
    M = D*(D+1)÷2
    N = D+M+1
    L ≤ N && error("The length of initial values should be larger than $(N).")
    F = zeros(L)
    X = ones(L, N)
    for i in 1:N-1
        p = ps[i]
        F[i] = f(p)
        j = 1
        for i1 in 1:D, i2 in i1:D
            X[i,j] = p[i1]*p[i2]
            j = j + 1
        end
        X[i,M+1:M+D] .= p
    end
    for _ in 1:n
        p = _recursion!(f, ps, F, X, L)
        push!(ps,p)
    end
    return ps
end

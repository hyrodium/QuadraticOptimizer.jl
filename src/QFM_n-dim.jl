function _recursion_qfm!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, F::Vector, X::Matrix, L::Integer) where {D}
    M = D*(D+1)÷2
    p = ps[end]
    i = mod(length(ps), 1:L)
    F[i] = fs[end]
    j = 1
    for i1 in 1:D, i2 in i1:D
        X[j,i] = p[i1]*p[i2]
        j = j + 1
    end
    X[M+1:M+D, i] .= p
    Y = pinv(X*X')*(X*F)
    a = Y[SOneTo(M)]
    b = SVector{D}(Y[M+1:M+D])
    c = Y[end]
    q = Quadratic(a,b,c)
    return center(q)
end

function optimize_qfm!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer) where D
    L = length(ps)
    M = D*(D+1)÷2
    N = D+M+1
    length(fs) == L ≥ N || error("The length of initial values should be larger than $(N).")
    F = zeros(L)
    X = ones(N, L)
    for i in 1:N-1
        p = ps[i]
        F[i] = fs[i]
        j = 1
        for i1 in 1:D, i2 in i1:D
            X[i,j] = p[i1]*p[i2]
            j = j + 1
        end
        X[M+1:M+D, i] .= p
    end
    for _ in 1:n
        p = _recursion_qfm!(f, ps, fs, F, X, L)
        push!(ps,p)
        push!(fs,f(p))
    end
    return ps, fs
end

function optimize_qfm(f, ps_init::Vector{<:SVector{D, <:Real}}, fs_init::Vector{<:Real}, n::Integer) where D
    ps = copy(ps_init)
    fs = copy(fs_init)
    return optimize_qfm!(f, ps, fs, n)
end

function optimize_qfm(f, ps_init::Vector{<:SVector{D, <:Real}}, n::Integer) where D
    ps = copy(ps_init)
    fs = f.(ps)
    return optimize_qfm!(f, ps, fs, n)
end

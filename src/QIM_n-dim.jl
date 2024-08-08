function _recursion_qim!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, F::MVector{N}, X::MMatrix{N,N}) where {D, N}
    i = mod(length(ps), 1:N)
    M = D*(D+1)÷2
    p = ps[end]
    F[i] = fs[end]
    j = 1
    for i1 in 1:D, i2 in i1:D
        X[j,i] = p[i1]*p[i2]
        j = j + 1
    end
    X[M+1:M+D, i] .= p
    Y = pinv(X') * F
    a = Y[SOneTo(M)]
    b = SVector{D}(Y[M+1:M+D])
    c = Y[end]
    q = Quadratic(a,b,c)
    return center(q)
end

function optimize_qim!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer) where D
    L = length(ps)
    M = D*(D+1)÷2
    N = D+M+1
    length(fs) == L == N || error("The length of initial values should be equal to $(N).")
    F = @MVector zeros(N)
    X = @MMatrix ones(N,N)
    for i in 1:N
        p = ps[i]
        F[i] = f(p)
        j = 1
        for i1 in 1:D, i2 in i1:D
            X[j,i] = p[i1]*p[i2]
            j = j + 1
        end
        X[M+1:M+D, i] .= p
    end
    for _ in 1:n
        p = _recursion_qim!(f, ps, fs, F, X)
        push!(ps,p)
        push!(fs,f(p))
    end
    return ps, fs
end

function optimize_qim(f, ps_init::Vector{<:SVector{D, <:Real}}, fs_init::Vector{<:Real}, n::Integer) where D
    ps = copy(ps_init)
    fs = copy(fs_init)
    return optimize_qim!(f, ps, fs, n)
end

function optimize_qim(f, ps_init::Vector{<:SVector{D, <:Real}}, n::Integer) where D
    ps = copy(ps_init)
    fs = f.(ps)
    return optimize_qim!(f, ps, fs, n)
end

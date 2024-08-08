function _recursion!(f, xs::Vector{<:Real}, F::MVector{3}, X::MVector{3})
    F[mod(length(xs), 1:3)] = f(xs[end])
    X[mod(length(xs), 1:3)] = xs[end]
    f₁ ,f₂, f₃ = F
    x₁, x₂, x₃ = X
    a₁ = f₁*(x₂-x₃)
    a₂ = f₂*(x₃-x₁)
    a₃ = f₃*(x₁-x₂)
    return (a₁*(x₂+x₃)+a₂*(x₃+x₁)+a₃*(x₁+x₂))/2(a₁+a₂+a₃)
end

function optimize_qim!(f, xs::Vector{<:Real}, n::Integer)
    length(xs) ≠ 3 && error("The length of initial values should be 3.")
    F = @MVector zeros(3)
    X = @MVector zeros(3)
    X[1] = xs[1]
    X[2] = xs[2]
    F[1] = f(xs[1])
    F[2] = f(xs[2])
    for _ in 1:n
        x = _recursion!(f, xs, F, X)
        push!(xs,x)
    end
    return xs
end

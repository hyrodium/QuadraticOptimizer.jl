# Developer Notes

## Symbols

| Symbol              | Type                        | Description                                | Remark            |
| :------------------ | :-------------------------- | :----------------------------------------- | :---------------- |
| `q`                 | `Quadratic{D,L,T}`          | Quadratic polynomial on ``\mathbb{R}^D``   |                   |
| `D`                 | `Integer`                   | Dimensions of the space ``\mathbb{R}^D``   |                   |
| `p`                 | `StaticVector{D}`           | Point in ``\mathbb{R}^D``                  |                   |
| `ps`                | `Vector{<:StaticVector{D}}` | Points in ``\mathbb{R}^D``                 |                   |
| `f`                 | `Function`, (`Any`)         | Objective on ``\mathbb{R}^D``              | Optimize = minimize in this package |
| `fs`                | `Vector{<:Real}`            | Objective values on `ps`                   | `f.(ps)`          |
| `ps_init`           | `Vector{<:StaticVector{D}}` | Initial points in ``\mathbb{R}^D``         | `ps[1:N]`         |
| `L`                 | `Integer`                   | Number of quadratic terms                  | `D*(D+1)/2`       |
| `M`, [`M(D)`](@ref) | `Integer`, `Function`       | Number of terms in quadratic polynomial    | `D+L+1`           |
| `N`                 | `Integer`                   | Number of input points                     | `length(ps_init)` |
| `n_iter`            | `Integer`                   | Number of maximum iterations               | Last arguments of [`optimize_qim`](@ref), [`optimize_qim!`](@ref), [`optimize_qfm`](@ref), and [`optimize_qfm!`](@ref) |
| `X`                 | `AbstractMatrix`            | Temporary matrix to store quadratic terms  | `size(X)==(M,N)` |
| `F`                 | `AbstractVector`            | Temporary vector to store objective values | `size(F)==(N)`   |

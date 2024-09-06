# Developer Notes

## Symbols

| Symbol      | Type                        | Description                                | Remark |
| :---------- | :-------------------------- | :----------------------------------------- | :----- |
| `q`,        | `Quadratic{D,L,T}`          | quadratic polynomial on ``\mathbb{R}^D``   |        |
| `D`         | `Integer`                   | Dimensions of the space ``\mathbb{R}^D``   |        |
| `p`         | `StaticVector`              | point in ``\mathbb{R}^D``                  |        |
| `ps`        | `Vector{<:StaticVector{D}}` | points in ``\mathbb{R}^D``                 |        |
| `f`         | `Function`, (`Any`)         | Objective on ``\mathbb{R}^D``              |        |
| `fs`        | `Vector{<:Real}`            | Objective values on `ps`                   |        |
| `L`         | `Integer`                   | Number of quadratic terms                  | `D*(D+1)÷2`       |
| `M`, `M(D)` | `Integer`, `Function`       | Number of terms of quadratic polynomial    | `D+L+1`           |
| `N`         | `Integer`                   | Number of input points                     | `length(ps_init)` |
| `n`         | `Integer`                   | Number of maximum iterations (argument)    |        |
| `X`         | `AbstractMatrix`            | Temporary matrix to store quadratic terms  |        |
| `F`         | `AbstractVector`            | Temporary vector to store objective values |        |

# QuadraticInterpolationMethod

[![Build Status](https://github.com/hyrodium/QuadraticInterpolationMethod.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hyrodium/QuadraticInterpolationMethod.jl/actions/workflows/CI.yml?query=branch%3Amain)


Quadratic interpolation method is an optimization method by interpolating given evaluation points with a quadratic polynomial.

```julia-repl
julia> using QuadraticInterpolationMethod: optimize!


julia> f(x) = sin(x) + x^2/10  # Function to minimize
f (generic function with 1 method)

julia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)
3-element Vector{Float64}:
  1.2
  0.1
 -2.2

julia> xs = copy(xs_init)  # Keep initial points
3-element Vector{Float64}:
  1.2
  0.1
 -2.2

julia> optimize!(f, xs, 10)  # Optimize 10 steps
13-element Vector{Float64}:
  1.2
  0.1
 -2.2
 -1.4980661244174434
 -1.2293686986818357
 -1.3061365335230135
 -1.3059492270208548
 -1.3064424808417185
 -1.3064400170208186
 -1.306440099017006
 -1.3066465256797584
 -1.306452471584103
 -1.3064400463690848
```

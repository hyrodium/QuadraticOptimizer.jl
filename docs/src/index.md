# QuadraticOptimizer.jl

A Julia implementation for quadratic interpolation method (QIM) and quadratic fitting method (QFM).

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hyrodium.github.io/QuadraticOptimizer.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hyrodium.github.io/QuadraticOptimizer.jl/dev)
[![Build Status](https://github.com/hyrodium/QuadraticOptimizer.jl/workflows/CI/badge.svg)](https://github.com/hyrodium/QuadraticOptimizer.jl/actions)
[![Coverage](https://codecov.io/gh/hyrodium/QuadraticOptimizer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/hyrodium/QuadraticOptimizer.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![QuadraticOptimizer Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/QuadraticOptimizer)](https://pkgs.genieframework.com?packages=QuadraticOptimizer)

## Quick start

```@repl quick_qim
using QuadraticOptimizer
f(x) = sin(x) + x^2/10  # Function to minimize
xs_init = [1.2, 0.1, -2.2]  # 3 initial points to construct a parabola
xs, fs = optimize_qim(f, xs_init, 10)  # Optimize 10 steps
```

```@example quick_qim
using Plots
pl = plot(f; xlims=(-5,5), color=:red3, label="objective")
plot!(pl, xs, fs; color=:blue3, label="iteration")
scatter!(pl, xs_init, f.(xs_init); color=:blue3, label="initial points")
```

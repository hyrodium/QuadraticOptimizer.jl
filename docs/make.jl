using Documenter
using QuadraticOptimizer
using QuadraticOptimizer: M
using Plots
using Random
using StaticArrays
using ForwardDiff

# Setup for doctests in docstrings
DocMeta.setdocmeta!(QuadraticOptimizer, :DocTestSetup, :(using QuadraticOptimizer, LinearAlgebra))

makedocs(;
    modules = [QuadraticOptimizer],
    format = Documenter.HTML(
        ansicolor=true,
        canonical = "https://hyrodium.github.io/QuadraticOptimizer.jl/stable/",
        assets = ["assets/favicon.ico", "assets/custom.css"],
        edit_link="main",
        repolink="https://github.com/hyrodium/QuadraticOptimizer.jl"
    ),
    pages = [
        "Home" => "index.md",
        "Developer Notes" => "dev-notes.md",
        "API" => "api.md",
    ],
    repo = "https://github.com/hyrodium/QuadraticOptimizer.jl/blob/{commit}{path}#L{line}",
    sitename = "QuadraticOptimizer.jl",
    authors = "hyrodium <hyrodium@gmail.com>",
)

deploydocs(
    repo = "github.com/hyrodium/QuadraticOptimizer.jl",
    push_preview = true,
    devbranch="main",
)

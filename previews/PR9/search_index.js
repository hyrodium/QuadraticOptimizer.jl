var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuadraticOptimizer]","category":"page"},{"location":"api/#QuadraticOptimizer.optimize_qfm!-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Vector{<:Real}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qfm!","text":"optimize_qfm!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Fitting Method (QFM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^d where f has been evaluated. This vector will be updated in-place during the optimization process.\nfs: A vector of function values corresponding to the points in ps. This vector will be updated in-place during the optimization process.\nn: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QFM, the last m points from ps and fs are used to fit a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:10]\n10-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n [0.34264792290134793, 0.2678704383989201]\n [0.515871408349502, 0.09002301604339691]\n [0.27265744579429385, 0.191562202596938]\n [0.4235912564725989, 0.4847023673932017]\n\njulia> ps = copy(ps_init);\n\njulia> fs = f.(ps);\n\njulia> optimize_qfm!(f, ps, fs, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qfm-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qfm","text":"optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, n::Integer)\n\nOptimize a function f using the Quadratic Fitting Method (QFM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^d where f has been evaluated.\nn: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QFM, the last m points from ps and fs are used to fit a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:10]\n10-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n [0.34264792290134793, 0.2678704383989201]\n [0.515871408349502, 0.09002301604339691]\n [0.27265744579429385, 0.191562202596938]\n [0.4235912564725989, 0.4847023673932017]\n\njulia> ps, fs = optimize_qfm(f, ps_init, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qfm-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Vector{<:Real}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qfm","text":"optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Fitting Method (QFM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^d where f has been evaluated.\nfs: A vector of function values corresponding to the points in ps.\nn: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QFM, the last m points from ps and fs are used to fit a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:10]\n10-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n [0.34264792290134793, 0.2678704383989201]\n [0.515871408349502, 0.09002301604339691]\n [0.27265744579429385, 0.191562202596938]\n [0.4235912564725989, 0.4847023673932017]\n\njulia> ps, fs = optimize_qfm(f, ps_init, f.(ps_init), 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim!-Tuple{Any, Vector{<:Real}, Vector{<:Real}, Integer}","page":"API","title":"QuadraticOptimizer.optimize_qim!","text":"optimize_qim!(f, xs::Vector{<:Real}, fs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nxs: A vector of points in mathbbR^d where f has been evaluated. This vector will be updated in-place during the optimization process.\nfs: A vector of function values corresponding to the points in xs. This vector will be updated in-place during the optimization process.\nn: The number of optimizing iterations. After execution, the length of xs will be m + n, where m = length(xs) before execution.\n\nnote: Note\nIn each step of the QIM, the last 3 points from xs and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending xs and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer\n\njulia> f(x) = sin(x) + x^2/10  # Function to minimize\nf (generic function with 1 method)\n\njulia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs = copy(xs_init);\n\njulia> fs = f.(xs);\n\njulia> optimize_qim!(f, xs, fs, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim!-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Vector{<:Real}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qim!","text":"optimize_qim!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^d where f has been evaluated. This vector will be updated in-place during the optimization process.\nfs: A vector of function values corresponding to the points in ps. This vector will be updated in-place during the optimization process.\nn: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QIM, the last m (==(D*(D+1)/2)) points from ps and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> ps = copy(ps_init);\n\njulia> fs = f.(ps);\n\njulia> optimize_qim!(f, ps, fs, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Tuple{Any, Vector{<:Real}, Integer}","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qfm(f, xs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nxs: A vector of points in mathbbR^d where f has been evaluated.\nn: The number of optimizing iterations. After execution, the length of xs will be m + n, where m = length(xs) before execution.\n\nnote: Note\nIn each step of the QIM, the last 3 points from xs and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending xs and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer\n\njulia> f(x) = sin(x) + x^2/10  # Function to minimize\nf (generic function with 1 method)\n\njulia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs = copy(xs_init)  # Keep initial points\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs, fs = optimize_qim(f, xs, 10);  # Optimize 10 steps\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Tuple{Any, Vector{<:Real}, Vector{<:Real}, Integer}","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qim(f, xs::Vector{<:Real}, fs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nxs: A vector of points in mathbbR^d where f has been evaluated.\nfs: A vector of function values corresponding to the points in xs.\nn: The number of optimizing iterations. After execution, the length of xs will be m + n, where m = length(xs) before execution.\n\nnote: Note\nIn each step of the QIM, the last 3 points from xs and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending xs and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer\n\njulia> f(x) = sin(x) + x^2/10  # Function to minimize\nf (generic function with 1 method)\n\njulia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs, fs = optimize_qim(f, xs_init, f.(xs_init), 10);  # Optimize 10 steps\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, n::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^d where f has been evaluated.\nn: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QIM, the last m (==(D*(D+1)/2)) points from ps and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> ps, fs = optimize_qfm(f, ps_init, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Vector{<:Real}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qim(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^d where f has been evaluated.\nfs: A vector of function values corresponding to the points in ps.\nn: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QIM, the last m (==(D*(D+1)/2)) points from ps and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> ps, fs = optimize_qim(f, ps_init, f.(ps_init), 20);\n\n\n\n\n\n","category":"method"},{"location":"#QuadraticOptimizer.jl","page":"Home","title":"QuadraticOptimizer.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia implementation for quadratic interpolation method (QIM) and quadratic fitting method (QFM).","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage) (Image: Aqua QA) (Image: QuadraticOptimizer Downloads)","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using QuadraticOptimizer\nf(x) = sin(x) + x^2/10  # Function to minimize\nxs_init = [1.2, 0.1, -2.2]  # 3 initial points to construct a parabola\nxs, fs = optimize_qim(f, xs_init, 10)  # Optimize 10 steps","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Plots\npl = plot(f; xlims=(-5,5), color=:red3, label=\"objective\")\nplot!(pl, xs, fs; color=:blue3, label=\"iteration\")\nscatter!(pl, xs_init, f.(xs_init); color=:blue3, label=\"initial points\")","category":"page"}]
}

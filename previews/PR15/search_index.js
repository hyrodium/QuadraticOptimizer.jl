var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [QuadraticOptimizer]","category":"page"},{"location":"api/#QuadraticOptimizer.Quadratic","page":"API","title":"QuadraticOptimizer.Quadratic","text":"Quadratic{D, L, T<:Real}\nQuadratic(a,b,c)\n\nA struct that represents quadratic polynomial.\n\nbeginaligned\nq(p) = frac12p^intercal A p + b^intercal p + c \nA = beginpmatrix\na_1      a_2       cdots  a_D \na_2      a_D+1   cdots  a_2D-1 \nvdots  vdots   ddots  vdots \na_D      a_2D-1  cdots  a_L\nendpmatrix\nendaligned\n\nExamples\n\njulia> q = Quadratic{2}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)\nQuadratic{2, 3, Float64}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)\n\njulia> q([1,2])\n7.7\n\n\n\n\n\n","category":"type"},{"location":"api/#QuadraticOptimizer.M-Tuple{Integer}","page":"API","title":"QuadraticOptimizer.M","text":"M(D::Integer) = (D+2)*(D+1)÷2\n\nCalculate the number of required points for QIM/QFM that is equal to the number of terms in D-dimensional quadratic polynomial.\n\nbeginaligned\nM(1) = 3 quad (textNumber of terms in  fraca_12 x^2 + b_1 x + c) \nM(2) = 6 quad (textNumber of terms in  fraca_12 x^2 + a_2 xy + fraca_32 y^2 + b_1 x + b_2 y + c) \nM(0) = 1 quad (textNumber of terms in  c)\nendaligned\n\nExamples\n\njulia> using QuadraticOptimizer: M\n\njulia> M(1)  # 3 points to interpolate a parabola\n3\n\njulia> M(2)  # 6 points to interpolate a 2-dimensional quadratic\n6\n\njulia> M(0)  # 1 point for constant, terms\n1\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.hessian-Tuple{Quadratic}","page":"API","title":"QuadraticOptimizer.hessian","text":"hessian(q::Quadratic)\n\nCalculate the (constant) hessian of the quadratic polynomial.\n\nExamples\n\njulia> using QuadraticOptimizer: hessian\n\njulia> using ForwardDiff\n\njulia> q = Quadratic{2}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)\nQuadratic{2, 3, Float64}([2.0, 1.0, 3.0], [-1.2, -2.3], 4.5)\n\njulia> hessian(q)\n2×2 StaticArrays.SHermitianCompact{2, Float64, 3} with indices SOneTo(2)×SOneTo(2):\n 2.0  1.0\n 1.0  3.0\n\njulia> ForwardDiff.hessian(q, [1,2])\n2×2 Matrix{Float64}:\n 2.0  1.0\n 1.0  3.0\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qfm!-Union{Tuple{T}, Tuple{D}, Tuple{Any, Vector{<:SVector{D, T}}, Vector{T}, Integer}} where {D, T<:Real}","page":"API","title":"QuadraticOptimizer.optimize_qfm!","text":"optimize_qfm!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Fitting Method (QFM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^D where f has been evaluated. This vector will be updated in-place during the optimization process.\nfs: A vector of function values corresponding to the points in ps. This vector will be updated in-place during the optimization process.\nn_iter: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QFM, the last m points from ps and fs are used to fit a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:10]\n10-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n [0.34264792290134793, 0.2678704383989201]\n [0.515871408349502, 0.09002301604339691]\n [0.27265744579429385, 0.191562202596938]\n [0.4235912564725989, 0.4847023673932017]\n\njulia> ps = copy(ps_init);\n\njulia> fs = f.(ps);\n\njulia> optimize_qfm!(f, ps, fs, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qfm-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qfm","text":"optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Fitting Method (QFM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^D where f has been evaluated.\nn_iter: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QFM, the last m points from ps and fs are used to fit a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:10]\n10-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n [0.34264792290134793, 0.2678704383989201]\n [0.515871408349502, 0.09002301604339691]\n [0.27265744579429385, 0.191562202596938]\n [0.4235912564725989, 0.4847023673932017]\n\njulia> ps, fs = optimize_qfm(f, ps_init, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qfm-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Vector{<:Real}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qfm","text":"optimize_qfm(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Fitting Method (QFM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^D where f has been evaluated.\nfs: A vector of function values corresponding to the points in ps.\nn_iter: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QFM, the last m points from ps and fs are used to fit a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:10]\n10-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n [0.34264792290134793, 0.2678704383989201]\n [0.515871408349502, 0.09002301604339691]\n [0.27265744579429385, 0.191562202596938]\n [0.4235912564725989, 0.4847023673932017]\n\njulia> ps, fs = optimize_qfm(f, ps_init, f.(ps_init), 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim!-Union{Tuple{T}, Tuple{Any, Vector{T}, Vector{T}, Integer}} where T<:Real","page":"API","title":"QuadraticOptimizer.optimize_qim!","text":"optimize_qim!(f, xs::Vector{<:Real}, fs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nxs: A vector of points in mathbbR^D where f has been evaluated. This vector will be updated in-place during the optimization process.\nfs: A vector of function values corresponding to the points in xs. This vector will be updated in-place during the optimization process.\nn_iter: The number of optimizing iterations. After execution, the length of xs will be m + n, where m = length(xs) before execution.\n\nnote: Note\nIn each step of the QIM, the last 3 points from xs and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending xs and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer\n\njulia> f(x) = sin(x) + x^2/10  # Function to minimize\nf (generic function with 1 method)\n\njulia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs = copy(xs_init);\n\njulia> fs = f.(xs);\n\njulia> optimize_qim!(f, xs, fs, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim!-Union{Tuple{T}, Tuple{D}, Tuple{Any, Vector{<:SVector{D, T}}, Vector{T}, Integer}} where {D, T<:Real}","page":"API","title":"QuadraticOptimizer.optimize_qim!","text":"optimize_qim!(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^D where f has been evaluated. This vector will be updated in-place during the optimization process.\nfs: A vector of function values corresponding to the points in ps. This vector will be updated in-place during the optimization process.\nn_iter: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QIM, the last M (==((D+2)*(D+1)/2)) points from ps and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> ps = copy(ps_init);\n\njulia> fs = f.(ps);\n\njulia> optimize_qim!(f, ps, fs, 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Tuple{Any, Vector{<:Real}, Integer}","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qim(f, xs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nxs: A vector of points in mathbbR^D where f has been evaluated.\nn_iter: The number of optimizing iterations. After execution, the length of xs will be m + n, where m = length(xs) before execution.\n\nnote: Note\nIn each step of the QIM, the last 3 points from xs and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending xs and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer\n\njulia> f(x) = sin(x) + x^2/10  # Function to minimize\nf (generic function with 1 method)\n\njulia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs = copy(xs_init)  # Keep initial points\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs, fs = optimize_qim(f, xs, 10);  # Optimize 10 steps\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Tuple{Any, Vector{<:Real}, Vector{<:Real}, Integer}","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qim(f, xs::Vector{<:Real}, fs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nxs: A vector of points in mathbbR^D where f has been evaluated.\nfs: A vector of function values corresponding to the points in xs.\nn_iter: The number of optimizing iterations. After execution, the length of xs will be m + n, where m = length(xs) before execution.\n\nnote: Note\nIn each step of the QIM, the last 3 points from xs and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending xs and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer\n\njulia> f(x) = sin(x) + x^2/10  # Function to minimize\nf (generic function with 1 method)\n\njulia> xs_init = [1.2, 0.1, -2.2]  # Initial points (3 points are required to construct parabola)\n3-element Vector{Float64}:\n  1.2\n  0.1\n -2.2\n\njulia> xs, fs = optimize_qim(f, xs_init, f.(xs_init), 10);  # Optimize 10 steps\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qim(f, ps::Vector{<:SVector{D, <:Real}}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^D where f has been evaluated.\nn_iter: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QIM, the last M (==((D+2)*(D+1)/2)) points from ps and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> ps, fs = optimize_qim(f, ps_init, 10);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.optimize_qim-Union{Tuple{D}, Tuple{Any, Vector{<:SVector{D, <:Real}}, Vector{<:Real}, Integer}} where D","page":"API","title":"QuadraticOptimizer.optimize_qim","text":"optimize_qim(f, ps::Vector{<:SVector{D, <:Real}}, fs::Vector{<:Real}, n_iter::Integer)\n\nOptimize a function f using the Quadratic Interpolation Method (QIM).\n\nArguments\n\nf: The objective function to be optimized.\nps: A vector of points in mathbbR^D where f has been evaluated.\nfs: A vector of function values corresponding to the points in ps.\nn_iter: The number of optimizing iterations. After execution, the length of ps will be m + n, where m = length(ps) before execution.\n\nnote: Note\nIn each step of the QIM, the last M (==((D+2)*(D+1)/2)) points from ps and fs are used to interpolate with a quadratic function. The method iteratively refines the points and function values, extending ps and fs with n additional points resulting from the optimization process.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5\nf (generic function with 1 method)\n\njulia> ps_init = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> ps, fs = optimize_qim(f, ps_init, f.(ps_init), 20);\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.quadratic_fitting-Union{Tuple{T}, Tuple{D}, Tuple{AbstractVector{<:StaticArray{Tuple{D}, T, 1}}, AbstractVector{T}}} where {D, T<:Real}","page":"API","title":"QuadraticOptimizer.quadratic_fitting","text":"quadratic_fitting(ps::AbstractVector{<:StaticVector{D, <:Real}}, fs::AbstractVector{<:Real}) where D\n\nCalculate quadratic fitting based on input values.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> q = Quadratic{2}([2,1,3], [1,2], 2)\nQuadratic{2, 3, Int64}([2, 1, 3], [1, 2], 2)\n\njulia> ps = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> quadratic_fitting(vcat(ps, ps), vcat(q.(ps).-1, q.(ps).+1))\nQuadratic{2, 3, Float64}([2.000000000387182, 0.9999999999890256, 2.999999999972893], [0.9999999997884512, 2.0000000000207985], 2.000000000052808)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuadraticOptimizer.quadratic_interpolation-Union{Tuple{T}, Tuple{D}, Tuple{AbstractVector{<:StaticArray{Tuple{D}, T, 1}}, AbstractVector{T}}} where {D, T<:Real}","page":"API","title":"QuadraticOptimizer.quadratic_interpolation","text":"quadratic_interpolation(ps::AbstractVector{<:StaticVector{D, <:Real}}, fs::AbstractVector{<:Real}) where D\n\nCalculate quadratic interpolation based on input values.\n\nExamples\n\njulia> using QuadraticOptimizer, StaticArrays, Random\n\njulia> Random.seed!(42);\n\njulia> q = Quadratic{2}([2,1,3], [1,2], 2)\nQuadratic{2, 3, Int64}([2, 1, 3], [1, 2], 2)\n\njulia> ps = [@SVector rand(2) for _ in 1:6]\n6-element Vector{SVector{2, Float64}}:\n [0.6293451231426089, 0.4503389405961936]\n [0.47740714343281776, 0.7031298490032014]\n [0.6733461456394962, 0.16589443479313404]\n [0.6134782250008441, 0.6683403279577278]\n [0.4570310908017041, 0.2993652953937611]\n [0.6611433726193705, 0.6394313620423493]\n\njulia> quadratic_interpolation(ps, q.(ps))\nQuadratic{2, 3, Float64}([1.9999999999997884, 0.999999999999996, 2.999999999999987], [1.0000000000001212, 2.0000000000000053], 1.999999999999966)\n\n\n\n\n\n","category":"method"},{"location":"#QuadraticOptimizer.jl","page":"Home","title":"QuadraticOptimizer.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Julia implementation for quadratic interpolation method (QIM) and quadratic fitting method (QFM).","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage) (Image: Aqua QA) (Image: QuadraticOptimizer Downloads)","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick start","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using QuadraticOptimizer\nf(x) = sin(x) + x^2/10  # Function to minimize\nxs_init = [1.2, 0.1, -2.2]  # 3 initial points to construct a parabola\nxs, fs = optimize_qim(f, xs_init, 10)  # Optimize 10 steps","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Plots\npl = plot(f; xlims=(-5,5), color=:red3, label=\"objective\")\nplot!(pl, xs, fs; color=:blue3, label=\"iteration\")\nscatter!(pl, xs_init, f.(xs_init); color=:blue3, label=\"initial points\")","category":"page"},{"location":"dev-notes/#Developer-Notes","page":"Developer Notes","title":"Developer Notes","text":"","category":"section"},{"location":"dev-notes/#Symbols","page":"Developer Notes","title":"Symbols","text":"","category":"section"},{"location":"dev-notes/","page":"Developer Notes","title":"Developer Notes","text":"Symbol Type Description Remark\nq, Quadratic{D,L,T} Quadratic polynomial on mathbbR^D \nD Integer Dimensions of the space mathbbR^D \np StaticVector{D} Point in mathbbR^D \nps Vector{<:StaticVector{D}} Points in mathbbR^D \nf Function, (Any) Objective on mathbbR^D Optimize = minimize in this package\nfs Vector{<:Real} Objective values on ps f.(ps)\nps_init Vector{<:StaticVector{D}} Initial points in mathbbR^D ps[1:N]\nL Integer Number of quadratic terms D*(D+1)/2\nM, M(D) Integer, Function Number of terms of quadratic polynomial D+L+1\nN Integer Number of input points length(ps_init)\nn_iter Integer Number of maximum iterations Last arguments of optimize_qim, optimize_qim!, optimize_qfm, and optimize_qfm!\nX AbstractMatrix Temporary matrix to store quadratic terms size(X)==(M,N)\nF AbstractVector Temporary vector to store objective values size(F)==(N)","category":"page"}]
}

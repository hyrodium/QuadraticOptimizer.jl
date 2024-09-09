using QuadraticOptimizer
using QuadraticOptimizer: center
using QuadraticOptimizer: hessian
using StaticArrays
using LinearAlgebra
using Random
using ForwardDiff
using Test
using Aqua

Aqua.test_all(QuadraticOptimizer)

@testset "Quadratic type" begin
    @testset "D = 1" begin
        D = 1
        L = D*(D+1)÷2
        a = SVector{L}(rand(L))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)
        A = hessian(q)
        T = Float64
        @test !iszero(q)
        @test iszero(q-q)
        @test -(2q+(+q)) ≈ 3(-q) == (-q)*3 == -3q == -2\q*6 == -6q/2
        @test iszero(Quadratic{D}(zeros(1), zeros(1), 0))
        @test iszero(Quadratic{D,T,L}(zeros(1), zeros(1), 0))
        @test isfinite(Quadratic{D,T,L}(zeros(1), zeros(1), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(1) ./ zeros(1), zeros(1), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(1), zeros(1) ./ zeros(1), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(1), zeros(1), 0 / 0))

        for _ in 1:10
            p = @SVector rand(D)
            @test q(center(q) + p) ≈ q(center(q) - p)
            @test (-q)(p) == -(q(p))
            @test (q+2q)(p) ≈ 3(q(p))
            @test q(p) ≈ p'*A*p/2 + b'*p + c
        end
    end

    @testset "D = 2" begin
        D = 2
        L = D*(D+1)÷2
        a = SVector{L}(rand(L))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)
        A = hessian(q)
        T = Float64
        @test !iszero(q)
        @test iszero(q-q)
        @test -(2q+(+q)) ≈ 3(-q) == (-q)*3 == -3q == -2\q*6 == -6q/2
        @test iszero(Quadratic{D}(zeros(3), zeros(2), 0))
        @test iszero(Quadratic{D,T,L}(zeros(3), zeros(2), 0))
        @test isfinite(Quadratic{D,T,L}(zeros(3), zeros(2), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(3) ./ zeros(3), zeros(2), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(3), zeros(2) ./ zeros(2), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(3), zeros(2), 0 / 0))

        for _ in 1:10
            p = @SVector rand(D)
            @test q(center(q) + p) ≈ q(center(q) - p)
            @test (-q)(p) == -(q(p))
            @test (q+2q)(p) ≈ 3(q(p))
            @test q(p) ≈ p'*A*p/2 + b'*p + c
        end
    end

    @testset "D = 3" begin
        D = 3
        L = D*(D+1)÷2
        a = SVector{L}(rand(L))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)
        A = hessian(q)
        T = Float64
        @test !iszero(q)
        @test iszero(q-q)
        @test -(2q+(+q)) ≈ 3(-q) == (-q)*3 == -3q == -2\q*6 == -6q/2
        @test iszero(Quadratic{D}(zeros(6), zeros(3), 0))
        @test iszero(Quadratic{D,T,L}(zeros(6), zeros(3), 0))
        @test isfinite(Quadratic{D,T,L}(zeros(6), zeros(3), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(6) ./ zeros(6), zeros(3), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(6), zeros(3) ./ zeros(3), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(6), zeros(3), 0 / 0))

        for _ in 1:10
            p = @SVector rand(D)
            @test q(center(q) + p) ≈ q(center(q) - p)
            @test (-q)(p) == -(q(p))
            @test (q+2q)(p) ≈ 3(q(p))
            @test q(p) ≈ p'*A*p/2 + b'*p + c
        end
    end
end

f1(p) = sin(p[1]) + p[1]^2/10
f2(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
q1 = Quadratic(SVector(1.2), SVector(4.2), -1.8)
q2 = Quadratic(SVector(1.2, -1.1/2, -0.8), SVector(4.2, -2.3), 1.7)
e1(p) = q1(SVector(p[1]))
e2(p) = q2(p)

@testset "exact" begin
    @testset "D = 1" begin
        xs_init = [1.2, 0.1, -2.2]
        ps_init = SVector{1}.(xs_init)
        xs, _ = optimize_qim(e1, xs_init, 1)
        ps, _ = optimize_qim(e1, ps_init, 1)
        @test ps[end] ≈ center(q1)
        @test xs[end] ≈ center(q1)[1]
    end

    @testset "D = 2" begin
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:6]
        ps, fs = optimize_qim(e2, ps_init, 1)
        @test ps[end] ≈ center(q2)

        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:10]
        ps, fs = optimize_qfm(e2, ps_init, 1)
        @test ps[end] ≈ center(q2)
    end
end

@testset "QIM" begin
    @testset "D = 1" begin
        xs_init = [1.2, 0.1, -2.2]
        ps_init = SVector{1}.(xs_init)
        xs, _ = optimize_qim(f1, xs_init, 6)
        ps, _ = optimize_qim(f1, ps_init, 6)
        @test all([p[1] for p in ps] .≈ xs)
        @test abs(ForwardDiff.derivative(f1, xs[end])) < 1e-6
        @test norm(ForwardDiff.gradient(f1, ps[end])) < 1e-6
        @test minimum(norm.(ForwardDiff.derivative.(f1, xs))) < 1e-6
        @test minimum(norm.(ForwardDiff.gradient.(f1, ps))) < 1e-5
    end

    @testset "D = 2" begin
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:6]
        ps, fs = optimize_qim(f2, ps_init, 20)
        @test length(ps) == 26
        @test norm(ForwardDiff.gradient(f2, ps[end])) < 1e-3
        @test norm(ForwardDiff.gradient(f2, ps[1])) > 1e-1
        @test minimum(norm.(ForwardDiff.gradient.(f2, ps))) < 1e-6
    end
end

@testset "QFM" begin
    @testset "D = 2" begin
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:10]
        ps, fs = optimize_qfm(f2, ps_init, 20)
        @test length(ps) == 30
        @test norm(ForwardDiff.gradient(f2, ps[end])) < 2e-3
        @test norm(ForwardDiff.gradient(f2, ps[1])) > 1e-1
        @test minimum(norm.(ForwardDiff.gradient.(f2, ps))) < 2e-4
    end
end

@testset "QIM-QFM" begin
    @testset "D = 2" begin
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:6]
        ps_qim, fs_qim = optimize_qim(f2, ps_init, 10)
        ps_qfm, fs_qfm = optimize_qfm(f2, ps_init, 10)
        @test ps_qim ≈ ps_qfm  atol=1e-5
        @test fs_qim ≈ fs_qfm
    end
end

@testset "interpolation, fitting" begin
    ps = [@SVector rand(2) for _ in 1:6]
    fs = e2.(ps)
    @test quadratic_interpolation(ps, fs) ≈ q2
    @test quadratic_fitting(ps, fs) ≈ q2

    ps = [@SVector rand(2) for _ in 1:10]
    fs = e2.(ps)
    @test_throws Exception quadratic_interpolation(ps, fs)
    @test quadratic_fitting(ps, fs) ≈ q2

    ps = [@SVector rand(2) for _ in 1:10000]
    fs = e2.(ps) + randn(10000)/1000
    @test_throws Exception quadratic_interpolation(ps, fs)
    @test !(quadratic_fitting(ps, fs) ≈ q2)
    @test norm((quadratic_fitting(ps, fs) - q2).a) < 1e-3
    @test norm((quadratic_fitting(ps, fs) - q2).b) < 1e-3
    @test norm((quadratic_fitting(ps, fs) - q2).c) < 1e-3
end

@testset "more precision types" begin
    @testset "Rational" begin
        f(p) = p[1]^3 - p[1]
        xs_init = big.([0//1, 1//2, 2//3])
        ps_init = SVector.(xs_init)
        xs, fxs = optimize_qim(f, xs_init, f.(xs_init), 10)
        ps, fps = optimize_qim(f, ps_init, f.(ps_init), 10)
        @test fxs isa Vector{Rational{BigInt}}
        @test fps isa Vector{Rational{BigInt}}
        @test fxs == fps
        @test SVector.(xs) == ps
        @test 1-3xs[end]^2 < 1e-30
    end

    @testset "BigFloat" begin
        f(p) = p[1]^3 - p[1]
        xs_init = big.([0/1, 1/2, 2/3])
        ps_init = SVector.(xs_init)
        xs, fxs = optimize_qim(f, xs_init, f.(xs_init), 10)
        ps, fps = optimize_qim(f, ps_init, f.(ps_init), 10)
        @test fxs isa Vector{BigFloat}
        @test fps isa Vector{BigFloat}
        @test fxs ≈ fps
        @test SVector.(xs) ≈ ps
        @test 1-3xs[end]^2 < 1e-30
    end
end

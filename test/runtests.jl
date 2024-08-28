using QuadraticOptimizer
using QuadraticOptimizer: center
using QuadraticOptimizer: interpolation
using QuadraticOptimizer: fitting
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
        M = D*(D+1)÷2
        a = SVector{M}(rand(M))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)

        for _ in 1:10
            x₁ = rand()
            @test q([x₁]) ≈ a[1]*x₁*x₁ + b[1]*x₁ + c
            @test q(center(q) + [x₁]) ≈ q(center(q) - [x₁])
        end
    end

    @testset "D = 2" begin
        D = 2
        M = D*(D+1)÷2
        a = SVector{M}(rand(M))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)

        for _ in 1:10
            x₁ = rand()
            x₂ = rand()
            @test q([x₁, x₂]) ≈ a[1]*x₁*x₁ + a[2]*x₁*x₂ + a[3]*x₂*x₂ + b[1]*x₁ + b[2]*x₂ + c
            @test q(center(q) + [x₁, x₂]) ≈ q(center(q) - [x₁, x₂])
        end
    end

    @testset "D = 3" begin
        D = 3
        M = D*(D+1)÷2
        a = SVector{M}(rand(M))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)

        for _ in 1:10
            x₁ = rand()
            x₂ = rand()
            x₃ = rand()
            @test q([x₁, x₂, x₃]) ≈
            (
                + a[1]*x₁*x₁ + a[2]*x₁*x₂ + a[3]*x₁*x₃
                             + a[4]*x₂*x₂ + a[5]*x₂*x₃
                                          + a[6]*x₃*x₃
                + b[1]*x₁
                + b[2]*x₂
                + b[3]*x₃
                + c
            )
            @test q(center(q) + [x₁, x₂, x₃]) ≈ q(center(q) - [x₁, x₂, x₃])
        end
    end
end

f1(p) = sin(p[1]) + p[1]^2/10
f2(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
q1 = Quadratic(SVector(1.2), SVector(4.2), -1.8)
q2 = Quadratic(SVector(1.2, -1.1, -0.8), SVector(4.2, -2.3), 1.7)
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
        ps, fs = optimize_qfm(f2, ps_init, 30)
        @test length(ps) == 40
        @test norm(ForwardDiff.gradient(f2, ps[end])) < 1e-3
        @test norm(ForwardDiff.gradient(f2, ps[1])) > 1e-1
        @test minimum(norm.(ForwardDiff.gradient.(f2, ps))) < 1e-4
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
    @test interpolation(ps, fs) ≈ q2
    @test fitting(ps, fs) ≈ q2

    ps = [@SVector rand(2) for _ in 1:10]
    fs = e2.(ps)
    @test_throws Exception interpolation(ps, fs)
    @test fitting(ps, fs) ≈ q2
end

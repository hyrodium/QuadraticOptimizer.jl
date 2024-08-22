using QuadraticOptimizer
using QuadraticOptimizer: center
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

        x₁ = rand()
        @test q([x₁]) ≈ a[1]*x₁*x₁ + b[1]*x₁ + c
    end

    @testset "D = 2" begin
        D = 2
        M = D*(D+1)÷2
        a = SVector{M}(rand(M))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)

        x₁ = rand()
        x₂ = rand()
        @test q([x₁, x₂]) ≈ a[1]*x₁*x₁ + a[2]*x₁*x₂ + a[3]*x₂*x₂ + b[1]*x₁ + b[2]*x₂ + c
    end

    @testset "D = 3" begin
        D = 3
        M = D*(D+1)÷2
        a = SVector{M}(rand(M))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)

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
    end
end

@testset "optimize_qim (exact)" begin
    @testset "D = 1" begin
        q = Quadratic(SVector(1.2), SVector(4.2), -1.8)
        f(p) = q(SVector(p[1]))
        xs_init = [1.2, 0.1, -2.2]
        ps_init = SVector{1}.(xs_init)
        xs, _ = optimize_qim(f, xs_init, 1)
        ps, _ = optimize_qim(f, ps_init, 1)
        @test ps[end] ≈ center(q)
        @test xs[end] ≈ center(q)[1]
    end

    @testset "D = 2" begin
        q = Quadratic(SVector(1.2, -1.1, -0.8), SVector(4.2, -2.3), 1.7)
        f(p) = q(p)
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:6]
        ps, fs = optimize_qim(f, ps_init, 1)
        @test ps[end] ≈ center(q)
    end
end

@testset "optimize_qim" begin
    @testset "D = 1" begin
        f(p) = sin(p[1]) + p[1]^2/10
        xs_init = [1.2, 0.1, -2.2]
        ps_init = SVector{1}.(xs_init)
        xs, _ = optimize_qim(f, xs_init, 6)
        ps, _ = optimize_qim(f, ps_init, 6)
        @test all([p[1] for p in ps] .≈ xs)
        @test abs(ForwardDiff.derivative(f, xs[end])) < 1e-6
        @test norm(ForwardDiff.gradient(f, ps[end])) < 1e-6
        @test minimum(norm.(ForwardDiff.derivative.(f, xs))) < 1e-6
        @test minimum(norm.(ForwardDiff.gradient.(f, ps))) < 1e-5
    end

    @testset "D = 2" begin
        f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:6]
        ps, fs = optimize_qim(f, ps_init, 20)
        @test length(ps) == 26
        @test norm(ForwardDiff.gradient(f, ps[end])) < 1e-3
        @test norm(ForwardDiff.gradient(f, ps[1])) > 1e-1
        @test minimum(norm.(ForwardDiff.gradient.(f, ps))) < 1e-6
    end
end

@testset "optimize_qfm (exact)" begin
    @testset "D = 2" begin
        q = Quadratic(SVector(1.2, -1.1, -0.8), SVector(4.2, -2.3), 1.7)
        f(p) = q(p)
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:10]
        ps, fs = optimize_qfm(f, ps_init, 1)
        @test ps[end] ≈ center(q)
    end
end

@testset "optimize_qfm" begin
    @testset "D = 2" begin
        f(p) = p[1]^2 + sin(p[1]) + 1.5p[2]^2 + sinh(p[2]) - p[1]*p[2]/5
        Random.seed!(42)
        ps_init = [@SVector rand(2) for _ in 1:10]
        ps, fs = optimize_qfm(f, ps_init, 30)
        @test length(ps) == 40
        @test norm(ForwardDiff.gradient(f, ps[end])) < 1e-3
        @test norm(ForwardDiff.gradient(f, ps[1])) > 1e-1
        @test minimum(norm.(ForwardDiff.gradient.(f, ps))) < 1e-4
    end
end

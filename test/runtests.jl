using QuadraticOptimizer
using QuadraticOptimizer: optimize!
using QuadraticOptimizer: center
using StaticArrays
using Test
using Aqua

Aqua.test_all(QuadraticOptimizer)

@testset "quadratic" begin
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

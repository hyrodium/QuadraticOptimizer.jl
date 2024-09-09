@testset "Quadratic type" begin
    @testset "constructor" begin
        # (a,b,c)
        ## Basic constructions
        @test Quadratic{1, Float32, 1}(SVector(1.0), SVector(2.0), 3.0) isa Quadratic{1,Float32,1}
        @test Quadratic{2, Float32, 3}(SVector(1.0, 2.0, 3.0), SVector(4.0, 5.0), 6.0) isa Quadratic{2,Float32,3}
        @test Quadratic{3, Float32, 6}(SVector(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), SVector(1.0, 2.0, 3.0), 4.0) isa Quadratic{3,Float32,6}
        @test Quadratic{1, Float32, 1}(SVector(1), SVector(2), 3) isa Quadratic{1,Float32,1}
        @test Quadratic{2, Float32, 3}(SVector(1, 2, 3), SVector(4, 5), 6) isa Quadratic{2,Float32,3}
        @test Quadratic{3, Float32, 6}(SVector(1, 2, 3, 4, 5, 6), SVector(1, 2, 3), 4) isa Quadratic{3,Float32,6}
        @test Quadratic{1, Float32}(SVector(1.0), SVector(2.0), 3.0) isa Quadratic{1,Float32,1}
        @test Quadratic{2, Float32}(SVector(1.0, 2.0, 3.0), SVector(4.0, 5.0), 6.0) isa Quadratic{2,Float32,3}
        @test Quadratic{3, Float32}(SVector(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), SVector(1.0, 2.0, 3.0), 4.0) isa Quadratic{3,Float32,6}
        @test Quadratic{1, Float32}(SVector(1), SVector(2), 3) isa Quadratic{1,Float32,1}
        @test Quadratic{2, Float32}(SVector(1, 2, 3), SVector(4, 5), 6) isa Quadratic{2,Float32,3}
        @test Quadratic{3, Float32}(SVector(1, 2, 3, 4, 5, 6), SVector(1, 2, 3), 4) isa Quadratic{3,Float32,6}
        @test Quadratic{1}([1.0], SVector(2.0), 3.0) isa Quadratic{1,Float64,1}
        @test Quadratic{2}([1.0, 2.0, 3.0], SVector(4.0, 5.0), 6.0) isa Quadratic{2,Float64,3}
        @test Quadratic{3}([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], SVector(1.0, 2.0, 3.0), 4.0) isa Quadratic{3,Float64,6}
        @test Quadratic{1}([1], SVector(2), 3) isa Quadratic{1,Int,1}
        @test Quadratic{2}([1, 2, 3], SVector(4, 5), 6) isa Quadratic{2,Int,3}
        @test Quadratic{3}([1, 2, 3, 4, 5, 6], SVector(1, 2, 3), 4) isa Quadratic{3,Int,6}
        @test Quadratic{1}([1.0], [2.0], 3.0) isa Quadratic{1,Float64,1}
        @test Quadratic{2}([1.0, 2.0, 3.0], [4.0, 5.0], 6.0) isa Quadratic{2,Float64,3}
        @test Quadratic{3}([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [1.0, 2.0, 3.0], 4.0) isa Quadratic{3,Float64,6}
        @test Quadratic{1}([1], [2], 3) isa Quadratic{1,Int,1}
        @test Quadratic{2}([1, 2, 3], [4, 5], 6) isa Quadratic{2,Int,3}
        @test Quadratic{3}([1, 2, 3, 4, 5, 6], [1, 2, 3], 4) isa Quadratic{3,Int,6}
        @test Quadratic(SVector(1.0), SVector(2.0), 3.0) isa Quadratic{1,Float64,1}
        @test Quadratic(SVector(1.0, 2.0, 3.0), SVector(4.0, 5.0), 6.0) isa Quadratic{2,Float64,3}
        @test Quadratic(SVector(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), SVector(1.0, 2.0, 3.0), 4.0) isa Quadratic{3,Float64,6}
        @test Quadratic(SVector(1), SVector(2), 3) isa Quadratic{1,Int,1}
        @test Quadratic(SVector(1, 2, 3), SVector(4, 5), 6) isa Quadratic{2,Int,3}
        @test Quadratic(SVector(1, 2, 3, 4, 5, 6), SVector(1, 2, 3), 4) isa Quadratic{3,Int,6}
        ## Method Error cases
        @test_throws MethodError Quadratic([1], [2], 3) isa Quadratic{1,Int,1}
        @test_throws MethodError Quadratic([1, 2, 3], [4, 5], 6) isa Quadratic{2,Int,3}
        @test_throws MethodError Quadratic([1, 2, 3, 4, 5, 6], [1, 2, 3], 4) isa Quadratic{3,Int,6}
        ## Input type promotion
        @test Quadratic(SVector(1,2,3), SVector(1//1, 1//2), 4f0) isa Quadratic{2,Float32}

        # (A,b,c)
        @test Quadratic(SMatrix{1,1}(1.0), SVector(2.0), 3.0) === Quadratic(SVector(1.0), SVector(2.0), 3.0)
        @test Quadratic(SMatrix{2,2}(1.0, 2.0, 2.0, 3.0), SVector(4.0, 5.0), 6.0) === Quadratic(SVector(1.0, 2.0, 3.0), SVector(4.0, 5.0), 6.0)
        @test Quadratic(SMatrix{3,3}(1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0), SVector(1.0, 2.0, 3.0), 4.0) === Quadratic(SVector(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), SVector(1.0, 2.0, 3.0), 4.0)

        # (c)
        @test Quadratic{2,Rational{Int},3}(1) isa Quadratic{2,Rational{Int},3}
        @test Quadratic{2,Rational{Int}}(1) isa Quadratic{2,Rational{Int},3}
        @test Quadratic{2}(1) isa Quadratic{2,Int,3}
        @test_throws MethodError Quadratic(1)
    end

    @testset "equality" begin
        @test Quadratic(SVector(1.0), SVector(2.0), 3.0) == Quadratic(SVector(1), SVector(2), 3)
        @test Quadratic(SVector(1.0, 2.0, 3.0), SVector(4.0, 5.0), 6.0) == Quadratic(SVector(1, 2, 3), SVector(4, 5), 6)
        @test Quadratic(SVector(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), SVector(1.0, 2.0, 3.0), 4.0) == Quadratic(SVector(1, 2, 3, 4, 5, 6), SVector(1, 2, 3), 4)
    end

    @testset "zero" begin
        # `Base.zero`
        @test zero(Quadratic{1}) === Quadratic(SVector(0.0), SVector(0.0), 0.0)
        @test zero(Quadratic{2}) === Quadratic(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0), 0.0)
        @test zero(Quadratic{3}) === Quadratic(SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0), 0.0)
        @test zero(Quadratic{1,Int}) === Quadratic(SVector(0), SVector(0), 0)
        @test zero(Quadratic{2,Int}) === Quadratic(SVector(0, 0, 0), SVector(0, 0), 0)
        @test zero(Quadratic{3,Int}) === Quadratic(SVector(0, 0, 0, 0, 0, 0), SVector(0, 0, 0), 0)
        @test_throws MethodError zero(Quadratic{1,Int,0})
        @test_throws MethodError zero(Quadratic{2,Int,0})
        @test_throws MethodError zero(Quadratic{3,Int,0})
        # `Base.iszero`
        @test iszero(zero(Quadratic{1}))
        @test iszero(zero(Quadratic{2}))
        @test iszero(zero(Quadratic{3}))
        @test iszero(zero(Quadratic{1,Int}))
        @test iszero(zero(Quadratic{2,Int}))
        @test iszero(zero(Quadratic{3,Int}))
        @test !iszero(Quadratic{1}(42))
        @test !iszero(Quadratic{2}(42))
        @test !iszero(Quadratic{3}(42))
        @test !iszero(Quadratic{1,Int}(42))
        @test !iszero(Quadratic{2,Int}(42))
        @test !iszero(Quadratic{3,Int}(42))
        # Values on zero-quadratic
        @test zero(Quadratic{1})(fill(42,1)) === 0.0
        @test zero(Quadratic{2})(fill(42,2)) === 0.0
        @test zero(Quadratic{3})(fill(42,3)) === 0.0
        @test zero(Quadratic{1,Int})(fill(42,1)) === 0.0
        @test zero(Quadratic{2,Int})(fill(42,2)) === 0.0
        @test zero(Quadratic{3,Int})(fill(42,3)) === 0.0
        @test zero(Quadratic{1,Rational{Int}})(fill(42,1)) === 0//1
        @test zero(Quadratic{2,Rational{Int}})(fill(42,2)) === 0//1
        @test zero(Quadratic{3,Rational{Int}})(fill(42,3)) === 0//1
        @test zero(Quadratic{1,Int})(fill(42//1,1)) === 0//1
        @test zero(Quadratic{2,Int})(fill(42//1,2)) === 0//1
        @test zero(Quadratic{3,Int})(fill(42//1,3)) === 0//1

        @testset "randomized (D = $D)" for D in 1:3
            for _ in 1:10
                L = D*(D+1)÷2
                a = SVector{L}(randn(L))
                b = SVector{D}(randn(D))
                c = randn()
                q = Quadratic(a,b,c)
                @test !iszero(q)
                @test iszero(q-q)
            end
        end
    end

    @testset "finite" begin
        for D in 1:3, T in (Float32, Float64)
            L = D*(D+1)÷2
            @test isfinite(Quadratic{D,T,L}(zeros(L), zeros(D), 0))
            @test !isfinite(Quadratic{D,T,L}(zeros(L) ./ zeros(L), zeros(D), 0))
            @test !isfinite(Quadratic{D,T,L}(zeros(L), zeros(D) ./ zeros(D), 0))
            @test !isfinite(Quadratic{D,T,L}(zeros(L), zeros(D), 0 / 0))
        end
    end

    @testset "D = $D" for D in 1:3
        L = D*(D+1)÷2
        a = SVector{L}(rand(L))
        b = SVector{D}(rand(D))
        c = rand()
        q = Quadratic(a,b,c)
        A = hessian(q)
        T = Float64
        @test -(2q+(+q)) ≈ 3(-q) == (-q)*3 == -3q == -2\q*6 == -6q/2

        for _ in 1:10
            p = @SVector rand(D)
            @test q(center(q) + p) ≈ q(center(q) - p)
            @test (-q)(p) == -(q(p))
            @test (q+2q)(p) ≈ 3(q(p))
            @test q(p) ≈ p'*A*p/2 + b'*p + c
        end
    end
end
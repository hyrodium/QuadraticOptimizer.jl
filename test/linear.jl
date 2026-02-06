@testset "Linear type" begin
    @testset "constructor" begin
        # (a, b)
        ## Basic constructions
        @test Linear{1, Float32}(SVector(1.0), 2.0) isa Linear{1,Float32}
        @test Linear{2, Float32}(SVector(1.0, 2.0), 3.0) isa Linear{2,Float32}
        @test Linear{3, Float32}(SVector(1.0, 2.0, 3.0), 4.0) isa Linear{3,Float32}
        @test Linear{1, Float32}(SVector(1), 2) isa Linear{1,Float32}
        @test Linear{2, Float32}(SVector(1, 2), 3) isa Linear{2,Float32}
        @test Linear{3, Float32}(SVector(1, 2, 3), 4) isa Linear{3,Float32}
        @test Linear{1}([1.0], 2.0) isa Linear{1,Float64}
        @test Linear{2}([1.0, 2.0], 3.0) isa Linear{2,Float64}
        @test Linear{3}([1.0, 2.0, 3.0], 4.0) isa Linear{3,Float64}
        @test Linear{1}([1], 2) isa Linear{1,Int}
        @test Linear{2}([1, 2], 3) isa Linear{2,Int}
        @test Linear{3}([1, 2, 3], 4) isa Linear{3,Int}
        @test Linear(SVector(1.0), 2.0) isa Linear{1,Float64}
        @test Linear(SVector(1.0, 2.0), 3.0) isa Linear{2,Float64}
        @test Linear(SVector(1.0, 2.0, 3.0), 4.0) isa Linear{3,Float64}
        @test Linear(SVector(1), 2) isa Linear{1,Int}
        @test Linear(SVector(1, 2), 3) isa Linear{2,Int}
        @test Linear(SVector(1, 2, 3), 4) isa Linear{3,Int}
        ## Method Error cases
        @test_throws MethodError Linear([1], 2)
        @test_throws MethodError Linear([1, 2], 3)
        @test_throws MethodError Linear([1, 2, 3], 4)
        ## Input type promotion
        @test Linear(SVector(1, 2), 3f0) isa Linear{2,Float32}
        @test Linear(SVector(1//1, 1//2), 4f0) isa Linear{2,Float32}

        # (b) constant only
        @test Linear{2,Rational{Int}}(1) isa Linear{2,Rational{Int}}
        @test Linear{2}(1) isa Linear{2,Int}
        @test Linear{3}(1.0) isa Linear{3,Float64}
        @test_throws MethodError Linear(1)
    end

    @testset "equality" begin
        @test Linear(SVector(1.0), 2.0) == Linear(SVector(1), 2)
        @test Linear(SVector(1.0, 2.0), 3.0) == Linear(SVector(1, 2), 3)
        @test Linear(SVector(1.0, 2.0, 3.0), 4.0) == Linear(SVector(1, 2, 3), 4)
    end

    @testset "zero" begin
        # `Base.zero`
        @test zero(Linear{1}) === Linear(SVector(0.0), 0.0)
        @test zero(Linear{2}) === Linear(SVector(0.0, 0.0), 0.0)
        @test zero(Linear{3}) === Linear(SVector(0.0, 0.0, 0.0), 0.0)
        @test zero(Linear{1,Int}) === Linear(SVector(0), 0)
        @test zero(Linear{2,Int}) === Linear(SVector(0, 0), 0)
        @test zero(Linear{3,Int}) === Linear(SVector(0, 0, 0), 0)
        @test zero(Linear{1}(2//1)) === Linear(SVector(0//1), 0//1)
        @test zero(Linear{2}(2//1)) === Linear(SVector(0//1, 0//1), 0//1)
        @test zero(Linear{3}(2//1)) === Linear(SVector(0//1, 0//1, 0//1), 0//1)
        @test zero(Linear{1,Int}(2//1)) === Linear(SVector(0), 0)
        @test zero(Linear{2,Int}(2//1)) === Linear(SVector(0, 0), 0)
        @test zero(Linear{3,Int}(2//1)) === Linear(SVector(0, 0, 0), 0)
        # `Base.iszero`
        @test iszero(zero(Linear{1}))
        @test iszero(zero(Linear{2}))
        @test iszero(zero(Linear{3}))
        @test iszero(zero(Linear{1,Int}))
        @test iszero(zero(Linear{2,Int}))
        @test iszero(zero(Linear{3,Int}))
        @test !iszero(Linear{1}(42))
        @test !iszero(Linear{2}(42))
        @test !iszero(Linear{3}(42))
        @test !iszero(Linear{1,Int}(42))
        @test !iszero(Linear{2,Int}(42))
        @test !iszero(Linear{3,Int}(42))
        # Values on zero-linear
        @test zero(Linear{1})(fill(42,1)) === 0.0
        @test zero(Linear{2})(fill(42,2)) === 0.0
        @test zero(Linear{3})(fill(42,3)) === 0.0
        @test zero(Linear{1,Int})(fill(42,1)) === 0
        @test zero(Linear{2,Int})(fill(42,2)) === 0
        @test zero(Linear{3,Int})(fill(42,3)) === 0
        @test zero(Linear{1,Rational{Int}})(fill(42,1)) === 0//1
        @test zero(Linear{2,Rational{Int}})(fill(42,2)) === 0//1
        @test zero(Linear{3,Rational{Int}})(fill(42,3)) === 0//1
        @test zero(Linear{1,Int})(fill(42//1,1)) === 0//1
        @test zero(Linear{2,Int})(fill(42//1,2)) === 0//1
        @test zero(Linear{3,Int})(fill(42//1,3)) === 0//1

        @testset "randomized test with (D = $D)" for D in 1:3
            for _ in 1:10
                a = SVector{D}(randn(D))
                b = randn()
                l = Linear(a, b)
                @test !iszero(l)
                @test iszero(l - l)
            end
        end
    end

    @testset "finite" begin
        for D in 1:3, T in (Float32, Float64)
            @test isfinite(Linear{D,T}(zeros(D), 0))
            @test !isfinite(Linear{D,T}(zeros(D) ./ zeros(D), 0))
            @test !isfinite(Linear{D,T}(zeros(D), 0 / 0))
        end
    end

    @testset "nan" begin
        for D in 1:3, T in (Float32, Float64)
            @test !isnan(Linear{D,T}(zeros(D), 0))
            @test isnan(Linear{D,T}(zeros(D) ./ zeros(D), 0))
            @test isnan(Linear{D,T}(zeros(D), 0 / 0))
        end
    end

    @testset "arithmetic" begin
        @test Linear{1}(4) !== Linear{2}(4)
        @test Linear{1}(4) != Linear{2}(4)
        @test Linear{1}(4) !== Linear{1}(4.0)
        @test Linear{1}(4) == Linear{1}(4.0)

        @test +Linear{2}(4) === Linear{2}(4)
        @test -Linear{1}(4) === Linear{1}(-4)
        @test 3 + Linear{4}(4) === Linear{4}(7)
        @test 3 - Linear{3}(4) === Linear{3}(-1)
        @test Linear{4}(4) + 3 === Linear{4}(7)
        @test Linear{3}(4) - 3 === Linear{3}(1)
        @test_throws DimensionMismatch Linear{2}(3) + Linear{4}(3)

        @testset "randomized test with (D = $D)" for D in 1:4
            for _ in 1:10
                a1 = SVector{D}(rand(D))
                a2 = SVector{D}(rand(D))
                b1 = rand()
                b2 = rand()
                l1 = Linear(a1, b1)
                l2 = Linear(a2, b2)
                @test l1 + l2 === l2 + l1 === Linear(a1+a2, b1+b2)
                @test l1 - l2 === -l2 + l1 === Linear(a1-a2, b1-b2)
                @test l1 + 1 === Linear(a1, b1+1)
                @test l1 + Linear{D}(1) === Linear(a1, b1+1)
                @test l1 + zero(l1) === l1
                @test zero(l1) + l1 === l1
                @test 2l1 === l1 + l1
                @test -l1 + l1 === l1 - l1 === zero(l1)

                for _ in 1:10
                    p = @SVector rand(D)
                    @test (-l1)(p) == -(l1(p))
                    @test (-l2)(p) == -(l2(p))
                    @test (l1 + 2l1 + l2/2)(p) ≈ 3(l1(p)) + l2(p)/2
                    @test l1(p) ≈ a1'*p + b1
                    @test l2(p) ≈ a2'*p + b2
                end
            end
        end
    end

    @testset "evaluation" begin
        l = Linear(SA[1.0, 2.0], 3.0)
        @test l([1.0, 1.0]) ≈ 6.0
        @test l([0.0, 0.0]) ≈ 3.0
        @test l([2.0, 3.0]) ≈ 1*2 + 2*3 + 3

        l2 = Linear(SA[1, 2, 3], 4)
        @test l2([1, 1, 1]) == 10
        @test l2([0, 0, 0]) == 4
    end
end

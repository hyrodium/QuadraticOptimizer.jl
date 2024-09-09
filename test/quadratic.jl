@testset "Quadratic type" begin
    @testset "D = $D" for D in 1:3
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
        @test iszero(Quadratic{D}(zeros(L), zeros(D), 0))
        @test iszero(Quadratic{D,T,L}(zeros(L), zeros(D), 0))
        @test isfinite(Quadratic{D,T,L}(zeros(L), zeros(D), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(L) ./ zeros(L), zeros(D), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(L), zeros(D) ./ zeros(D), 0))
        @test !isfinite(Quadratic{D,T,L}(zeros(L), zeros(D), 0 / 0))

        for _ in 1:10
            p = @SVector rand(D)
            @test q(center(q) + p) ≈ q(center(q) - p)
            @test (-q)(p) == -(q(p))
            @test (q+2q)(p) ≈ 3(q(p))
            @test q(p) ≈ p'*A*p/2 + b'*p + c
        end
    end
end
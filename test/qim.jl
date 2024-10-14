@testset "QIM" begin
    @testset "D = 1 (ps-xs)" begin
        xs_init = [1.2, 0.1, -2.2]
        ps_init = SVector{1}.(xs_init)
        xs, _ = optimize_qim(F, xs_init, 5)
        ps, _ = optimize_qim(F, ps_init, 5)
        @test all([p[1] for p in ps] .≈ xs)
    end

    @testset "D = $(D)" for D in 1:5
        Random.seed!(1)
        ps_init = [@SVector rand(D) for _ in 1:M(D)]
        ps, fs = optimize_qim(F, ps_init, M(D))
        @test ps_init == ps[1:M(D)]
        @test length(ps) == 2M(D)
        @test norm(ForwardDiff.gradient(F, ps[end])) < 1e-2
        @test norm(ForwardDiff.gradient(F, ps[1])) > 3e-1
        @test minimum(norm.(ForwardDiff.gradient.(F, ps))) < 3e-3
    end

    @testset "exact" begin
        @testset "D = 1 (ps-xs)" begin
            D = 1
            xs_init = [1.2, 0.1, -2.2]
            ps_init = SVector{1}.(xs_init)
            xs, _ = optimize_qim(x->Q(Val(D))(SVector(x)), xs_init, 1)
            ps, _ = optimize_qim(Q(Val(D)), ps_init, 1)
            @test ps[end] ≈ center(Q(Val(D)))
            @test xs[end] ≈ center(Q(Val(D)))[1]
        end
    
        @testset "D = $(D)" for D in 1:5
            Random.seed!(42)
            ps_init = [@SVector rand(D) for _ in 1:M(D)]
            ps, fs = optimize_qim(Q(Val(D)), ps_init, 1)
            @test ps[end] ≈ center(Q(Val(D)))
        end
    end

    @testset "interpolation" begin
        q2 = Quadratic(SVector(1.2, -1.1/2, -0.8), SVector(4.2, -2.3), 1.7)
        ps = [@SVector rand(2) for _ in 1:6]
        fs = q2.(ps)
        @test quadratic_interpolation(ps, fs) ≈ q2

        ps = [@SVector rand(2) for _ in 1:10]
        fs = q2.(ps)
        @test_throws Exception quadratic_interpolation(ps, fs)

        ps = [@SVector rand(2) for _ in 1:10000]
        fs = q2.(ps) + randn(10000)/1000
        @test_throws Exception quadratic_interpolation(ps, fs)
    end
end

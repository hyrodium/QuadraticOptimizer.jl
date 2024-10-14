@testset "QIM-QFM" begin
    @testset "D = $(D)" for D in 1:5
        Random.seed!(42)
        ps_init = [@SVector rand(D) for _ in 1:M(D)]
        ps_qim, fs_qim = optimize_qim(F, ps_init, 4)
        ps_qfm, fs_qfm = optimize_qfm(F, ps_init, 4)
        @test ps_qim ≈ ps_qfm  atol=1e-5
        @test fs_qim ≈ fs_qfm
    end

    @testset "more precision types" begin
        @testset "Rational" begin
            xs_init = big.([0//1, 1//2, 2//3])
            ps_init = SVector.(xs_init)
            xs_qim, fxs_qim = optimize_qim(G, xs_init, G.(xs_init), 12)
            ps_qim, fps_qim = optimize_qim(G, ps_init, G.(ps_init), 12)
            xs_qfm, fxs_qfm = optimize_qfm(G, xs_init, G.(xs_init), 12)
            ps_qfm, fps_qfm = optimize_qfm(G, ps_init, G.(ps_init), 12)
            @test fxs_qim isa Vector{Rational{BigInt}}
            @test fps_qim isa Vector{Rational{BigInt}}
            @test fxs_qfm isa Vector{Rational{BigInt}}
            @test fps_qfm isa Vector{Rational{BigInt}}
            @test fxs_qim == fps_qim
            @test_broken fxs_qim == fxs_qfm
            @test fxs_qim == fps_qfm
            @test SVector.(xs_qim) == ps_qim
            @test ForwardDiff.derivative(G, xs_qim[end]) < 1e-20
        end

        @testset "BigFloat" begin
            xs_init = BigFloat.([0//1, 1//2, 2//3])
            ps_init = SVector.(xs_init)
            xs_qim, fxs_qim = optimize_qim(G, xs_init, G.(xs_init), 12)
            ps_qim, fps_qim = optimize_qim(G, ps_init, G.(ps_init), 12)
            xs_qfm, fxs_qfm = optimize_qfm(G, xs_init, G.(xs_init), 12)
            ps_qfm, fps_qfm = optimize_qfm(G, ps_init, G.(ps_init), 12)
            @test fxs_qim isa Vector{BigFloat}
            @test fps_qim isa Vector{BigFloat}
            @test fxs_qfm isa Vector{BigFloat}
            @test fps_qfm isa Vector{BigFloat}
            @test fxs_qim ≈ fps_qim
            @test_broken fxs_qim ≈ fxs_qfm
            @test fxs_qim ≈ fps_qfm
            @test SVector.(xs_qim) ≈ ps_qim
            @test ForwardDiff.derivative(G, xs_qim[end]) < 1e-20
        end
    end
end

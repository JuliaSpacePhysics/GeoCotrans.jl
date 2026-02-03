@testitem "trace" begin
    using OrdinaryDiffEqTsit5
    using Dates
    using LinearAlgebra
    t = DateTime(2020, 1, 1)

    @testset "Tracing in both directions" begin
        pos = [3.0, 0.0, 0.5]
        in = (GEO(), Cartesian3())

        # Trace parallel to B
        sol_fwd = trace(pos, t, Tsit5(); dir = 1, in)

        # Trace anti-parallel to B
        sol_bwd = trace(pos, t, Tsit5(); dir = -1, in)

        # Both should reach inner boundary for dipole-like field
        r_fwd_end = norm(sol_fwd.u[end])
        r_bwd_end = norm(sol_bwd.u[end])
        @test r_fwd_end ≈ 1.0 rtol = 0.1
        @test r_bwd_end ≈ 1.0 rtol = 0.1

        # The endpoints should be in different hemispheres (approximately)
        z_fwd = sol_fwd.u[end][3]
        z_bwd = sol_bwd.u[end][3]
        @test sign(z_fwd) != sign(z_bwd) || abs(z_fwd) < 0.5 || abs(z_bwd) < 0.5
    end

    @testset "Outer boundary reached" begin
        pos = GEO(-5.0, 0.0, 3.0)
        sol = trace(pos, t, Tsit5(); rlim = 6.0, dir = -1)

        # Should eventually hit outer boundary
        @test length(sol.u) > 1
        r_end = norm(sol.u[end])
        @test 5.9 < r_end <= 6.1
    end

    @testset "CoordinateVector input interface" begin
        pos = GEO(3.0, 0.0, 0.0)
        sol = trace(pos, t, Tsit5(); rlim = 10.0, r0 = 1.0)

        @test norm(sol.u[end]) ≈ 1.0 rtol = 0.1
        @test length(sol.u) > 10
    end
end


@testitem "TsyganenkoModels field line tracing" begin
    using Dates
    using TsyganenkoModels, GeoCotrans
    using OrdinaryDiffEqTsit5
    using Chairmarks

    # Test case from https://github.com/tsssss/geopack/blob/master/geopack/test_geopack1.md
    # Note: their TS04 case is not up to date
    time = DateTime(2001, 1, 1, 2, 3, 4)
    r_gsm = GSM(-5.1, 0.3, 2.8)
    param = (; pdyn = 2.0, dst = -87.0, byimf = 2.0, bzimf = -5.0)

    # Expected footpoints from geopack test (GSM coordinates)
    footpoints = (;
        T89 = GSM(-0.7218581, 0.032782, 0.8293646),
        T96 = GSM(-0.7289211, 0.0353549, 0.8230575),
        T01 = GSM(-0.7264082, 0.0389479, 0.8251144),
        # TS04 = GSM(-0.7153766644654653, 0.024213930362513274, 0.8352541020689167),
    )

    models = (;
        T89 = TsyIGRF(T89(iopt = 2)),
        T96 = TsyIGRF(T96(param)),
        T01 = TsyIGRF(T01(param)),
        # TS04 = TsyIGRF(TS04(param)),
    )
    benchmark = @b trace($r_gsm, $time, Tsit5(); model = TsyIGRF(), r0 = 1.1)
    if VERSION >= v"1.12"
        @test benchmark.allocs < 150
    end
    @info benchmark

    for name in keys(models)
        model = models[name]
        sol = trace(r_gsm, time, Tsit5(); model, r0 = 1.1, rlim = 10.0)
        footpoint = sol.u[end]
        @test footpoint ≈ footpoints[name] rtol = 1.0e-3
    end
end

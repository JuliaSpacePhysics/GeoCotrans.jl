using TestItems, TestItemRunner
@run_package_tests

# https://github.com/spedas/pyspedas/blob/master/pyspedas/cotrans_tools/tests/cotrans.py

@testitem "CoordinateVector" begin
    using Dates
    using GeoCotrans: getcsys, GDZ, GEO, Geodetic

    𝐫 = GDZ(0, 60, 5)
    @test getcsys(𝐫) == getcsys(GDZ) == (GEO(), Geodetic())
    @test 𝐫 .* 2.0 isa CoordinateVector{GEO, Geodetic}

    t = Date(2021, 3, 28)
    𝐫_t = GDZ(0, 60, 5, t)
    @test getcsys(𝐫_t) == getcsys(GDZ)
    @test (𝐫_t .* 2.0).t == nothing
end

@testitem "gse2gsm" begin
    using Dates
    using DimensionalData

    t = Date(2021, 1, 1)
    gse = GSE(1, 2, 3)
    gse_t = GSE(1, 2, 3, t)
    gsm = gse2gsm(gse, t)
    @test gsm isa CoordinateVector{GSM}
    @test gsm.t == t
    @test GSM(gse, t) === gsm
    @test GSM(gse_t) === gsm
    @test_throws AssertionError gse2gsm(GSM(1, 2, 3), t)

    # Test GSE->GSM transformation
    # Data from Python test case
    data = [
        [775.0, 10.0, -10.0],
        [245.0, -102.0, 251.0],
        [121.0, 545.0, -1.0],
        [304.65, -205.3, 856.1],
        [464.34, -561.55, -356.22],
    ]

    # Convert Unix timestamps to DateTime objects
    timestamps = [
        unix2datetime(1577308800),
        unix2datetime(1577112800),
        unix2datetime(1577598800),
        unix2datetime(1577608800),
        unix2datetime(1577998800),
    ]

    # Create a DimArray with the test data
    da = DimArray(permutedims(hcat(data...)), (Ti(timestamps), Y(1:3)))
    gsm_da = gse2gsm(da)

    # Check that the transformation worked correctly
    expected_gsm = [775.0, 11.70357713, -7.93890939]
    @test gsm_da[Ti = 1] ≈ expected_gsm rtol = 1.0e-6
end

@testitem "Geodetic to Cartesian" begin
    using GeoCotrans: gdz2car, gdz2sph, car2sph
    using LinearAlgebra
    lat, lon, alt = 60, 5, 0
    𝐫 = gdz2car(lat, lon, alt)
    @test 𝐫 ≈ [0.4998961951811961, 0.043735250018097305, 0.8633345576875061]
    @test norm(𝐫) ≈ 1.0 rtol = 1.0e-2
    @test gdz2sph(lat, lon, alt) ≈ car2sph(𝐫)
end


@testitem "Validation with CoordRefSystems" begin
    using Dates
    using GeoCotrans: R🜨
    using CoordRefSystems
    import CoordRefSystems as CRS
    using Unitful

    lat, lon, alt = 60, 5, 0
    c1 = GeocentricLatLonAlt{WGS84{2139}}(lat, lon, alt)
    c2 = convert(CRS.Spherical, convert(Cartesian, c1))

    rθϕ = gdz2sph(lat, lon, alt)
    @test c2.r ≈ rθϕ[1] * R🜨 * Unitful.km rtol = 1.0e-5
    @test_broken c2.r ≈ rθϕ[1] * R🜨 * Unitful.km rtol = 1.0e-6
    @test c2.θ ≈ rθϕ[2] rtol = 1.0e-2 # our implementation is not as accurate
    @test_broken c2.θ ≈ rθϕ[2] rtol = 1.0e-5
    @test c2.ϕ ≈ rθϕ[3]
end

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(GeoCotrans)
end

@testitem "JET" begin
    using JET
    @test_call GeoCotrans.workload()
    @test_opt GeoCotrans.workload()
end

@testitem "trace_field_line" begin
    using Dates
    using LinearAlgebra
    using OrdinaryDiffEq

    t = DateTime(2020, 1, 1)

    # Basic tracing from dayside equator (GEO coordinates)
    pos = [4.0, 0.0, 0.0]
    sol = trace_field_line(pos, t, Tsit5(); rlim=10.0, r0=1.0)

    # Should have multiple points
    @test length(sol.u) > 10
    # First point should be the starting position
    @test sol.u[1] ≈ pos
    # Last point should be at r ≈ r0 (inner boundary)
    r_end = sqrt(sum(sol.u[end].^2))
    @test r_end ≈ 1.0 rtol=0.1

    # Test both directions
    sol_fwd = trace_field_line([3.0, 0.0, 0.5], t, Tsit5(); dir=1, rlim=10.0, r0=1.0)
    sol_bwd = trace_field_line([3.0, 0.0, 0.5], t, Tsit5(); dir=-1, rlim=10.0, r0=1.0)
    @test length(sol_fwd.u) > 5
    @test length(sol_bwd.u) > 5

    # Test CoordinateVector input
    pos_geo = GEO(3.0, 0.0, 0.0)
    sol_cv = trace_field_line(pos_geo, t, Tsit5(); rlim=10.0, r0=1.0)
    @test length(sol_cv.u) > 5

    # Test FieldLineProblem directly
    prob = FieldLineProblem([3.0, 0.0, 0.0], (0.0, 50.0), t)
    cb = FieldLineCallback(r0=1.0, rlim=10.0)
    sol_prob = solve(prob, Tsit5(); callback=cb)
    @test length(sol_prob.u) > 5
end

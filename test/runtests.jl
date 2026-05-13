using Test
using TestItems, TestItemRunner
@run_package_tests

# https://github.com/spedas/pyspedas/blob/master/pyspedas/cotrans_tools/tests/test_cotrans.py
# https://github.com/tsssss/geopack/blob/master/geopack/test_geopack1.py

const RUN_JET_TESTS = isempty(VERSION.prerelease)

if RUN_JET_TESTS
    using Pkg; Pkg.add("JET"); Pkg.instantiate()
    @testitem "JET static analysis" begin
        using JET
        @test_call GeoCotrans.workload()
        @test_opt GeoCotrans.workload()
    end
end

@testitem "Frames and Representations" begin
    using GeoCotrans: getcsys, representation
    # Default coordinate representation is Cartesian3
    @test representation(GEO()) == Cartesian3()
end

@testitem "CoordinateVector" begin
    using Dates
    using GeoCotrans: getcsys, GDZ, GEO, Geodetic, representation

    𝐫 = GDZ(0, 60, 5)
    @test getcsys(𝐫) == getcsys(GDZ) == (GEO(), Geodetic())
    @test 𝐫 .* 2.0 isa GeoCotrans.CoordinateVector{GEO, Geodetic}
    @test representation(𝐫) == Geodetic()
end

@testitem "gse2gsm" begin
    using Dates
    using DimensionalData

    t = Date(2021, 1, 1)
    gse = [1, 2, 3]
    gsm = gse2gsm(gse, t)
    @test gsm isa AbstractVector
    @test size(gsm) == (3,)
    @test_throws AssertionError gse2gsm(GeoCotrans.CoordinateVector{GSM, Cartesian3}(1, 2, 3), t)

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

    @testset "model evaluation for DimArray" begin
        @test_nowarn igrf(gsm_da; in = (GSM(), Cartesian3()), scale = 1 / 6371.2)
    end
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

@testitem "geo2mag" begin
    using Dates
    using GeoCotrans: calc_dipole_geo, calc_dipole_angle, get_igrf_coeffs

    t = DateTime(2001, 1, 1, 2, 3, 4)

    # Test that dipole direction in MAG is along Z-axis
    dipole_geo = calc_dipole_geo(t)
    dipole_mag = rotation(GEO, MAG, t) * dipole_geo
    @test dipole_mag ≈ [0, 0, 1]

    # Test CoordinateVector transformation
    geo = [2.8011944117533565, -2.4048913761357267, 4.5066403602406275]
    mag = geo2mag(geo, t)
    r_mag = [2.298686529948157, 1.8997853069109853, 5.004683409035982]
    @test mag ≈ r_mag
    @test mag2geo(mag, t) ≈ geo
end

@testitem "mlt" begin
    using Dates
    using GeoCotrans: get_mlt
    using DimensionalData
    using Chairmarks

    r = [2.195517156287977, 2.834061428571752, 0.34759070278576953]
    t = DateTime("2015-02-02T06:12:43")
    IRBEM_MLT = 9.56999052595853
    @test get_mlt(r, t) ≈ IRBEM_MLT rtol = 1.0e-4

    times = t .+ Minute.(0:2)
    A = stack([r, r, r]; dims = 1)
    da = DimArray(A, (Ti(times), Y(1:3)))
    @test get_mlt(da) ≈ get_mlt(A, times; dim = 1)

    @test (@b get_mlt($r, $t)).allocs == 0
end

@testitem "gsm2sm" begin
    using Dates
    using GeoCotrans: dipole_tilt, calc_dipole_gei

    t = DateTime(2001, 1, 1, 2, 3, 4)
    M = rotation(GSM, SM, t)
    # Test that dipole in GSM transforms to Z-axis in SM
    dipole_gsm = rotation(GEI, GSM, t) * calc_dipole_gei(t)
    dipole_sm = M * dipole_gsm
    @test dipole_sm ≈ [0, 0, 1]
    # Verify dipole tilt angle is consistent with dipole in GSM
    μ = dipole_tilt(t)
    @test dipole_gsm ≈ [sin(μ), 0, cos(μ)]

    # Test CoordinateVector transformation
    gsm = [-5.1, 0.3, 2.8]
    sm = gsm2sm(gsm, t)
    sm_geopack_true = [-2.9670092644479498, 0.3, 5.004683409035982]
    @test sm ≈ sm_geopack_true rtol = 1.0e-7
    @test sm2gsm(sm, t) ≈ gsm
    # Test roundtrip transformation
    @test gsm2sm(sm2gsm(sm, t), t) ≈ sm
end

# The most expensive transformation
@testitem "mag2sm and mag2gei" begin
    using GeoCotrans: mag2sm
    using Chairmarks
    using Dates

    t = DateTime(2001, 1, 1, 2, 3, 4)
    mag = [1.0, 2.0, 3.0]
    sm = mag2sm(mag, t)

    @info @b mag2sm($mag, $t)
    @test sm2mag(sm, t) ≈ mag
    @test gei2mag(mag2gei(mag, t), t) ≈ mag
end

@testitem "geo2sm" begin
    using Dates

    t = DateTime(2001, 1, 1, 2, 3, 4)
    geo = [1.0, 2.0, 3.0]
    @test sm2geo(geo2sm(geo, t), t) ≈ geo
end

@testitem "rotation graph" begin
    using Dates
    using LinearAlgebra

    t = DateTime(2001, 1, 1, 2, 3, 4)

    # identity
    @test rotation(GEI, GEI, t) == I
    @test rotation(GSM, GSM, t) == I

    # inverse via transpose round-trips to identity (direct edges)
    for (F1, F2) in ((GEI, GSM), (GSE, GSM), (GEO, MAG), (GSM, SM), (GEI, GEO))
        R = rotation(F1, F2, t)
        @test R * rotation(F2, F1, t) ≈ I
        @test rotation(F2, F1, t) ≈ transpose(R)
    end

    # all chain pairs round-trip to identity
    for (F1, F2) in ((GEI, SM), (GEO, GSM), (GEO, SM), (GEI, MAG), (SM, MAG))
        @test rotation(F1, F2, t) * rotation(F2, F1, t) ≈ I
    end
end

@testitem "transform interface" begin
    using Dates
    using LinearAlgebra

    t = DateTime(2001, 1, 1, 2, 3, 4)
    R = rotation(GEI, GSM, t)

    # raw vector with explicit (to, from)
    x = [1.0, 2.0, 3.0]
    @test transform(GSM, GEI, x, t) ≈ R * x

    # CoordinateVector: frame inferred via internal vector type
    gei = GeoCotrans.CoordinateVector{GEI, Cartesian3}(1.0, 2.0, 3.0)
    gsm = transform(GSM, gei, t)
    @test gsm isa GeoCotrans.CoordinateVector{GSM, Cartesian3}
    @test gsm ≈ gei2gsm(gei, t)

    @test transform(GSM, gei, t) ≈ gsm

    # frame mismatch is caught
    @test_throws AssertionError transform(GSM, GSE, gei, t)

    # identity
    @test transform(GEI, GEI, x, t) ≈ x
    @test transform(GEI, gei, t) ≈ gei

    # matrix path: scalar time builds rotation once
    A = rand(3, 5)
    @test transform(GSM, GEI, A, t) ≈ R * A

    # matrix path: per-sample times
    ts = [t + Hour(i) for i in 0:4]
    B = transform(GSM, GEI, A, ts; dims = 2)
    @test B ≈ gei2gsm(A, ts; dims = 2)

    # matrix path: dims = 1
    A1 = permutedims(A)
    @test transform(GSM, GEI, A1, ts; dims = 1) ≈ permutedims(B)

    # Pair syntax: from=>to
    @test transform(GEI => GSM, x, t) ≈ transform(GSM, GEI, x, t)
    @test transform(GEI => GSM, gei, t) ≈ gsm
    @test transform(GEI => GSM, A, t) ≈ R * A
    @test transform(GEI => GSM, A, ts; dims = 2) ≈ B
end

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(GeoCotrans)
end

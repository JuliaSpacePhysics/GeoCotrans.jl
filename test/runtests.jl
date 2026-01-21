using TestItems, TestItemRunner
@run_package_tests

# https://github.com/spedas/pyspedas/blob/master/pyspedas/cotrans_tools/tests/cotrans.py

@testitem "CoordinateVector" begin
    𝐫 = GDZ(0, 60, 5)
    @test getcsys(𝐫) == GDZ()
    @test 𝐫 .* 2.0 isa CoordinateVector{GDZ}
end

@testitem "gse2gsm" begin
    using Dates
    using DimensionalData

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
    @test isapprox(gsm_da[Ti = 1], expected_gsm)
end

@testitem "Validation with AstroLib" begin
    using Dates
    for yr in 2010:2025
        dt = Date(yr, 1, 1)
        # test GST
        @test isapprox(GeoCotrans.calculate_gst(dt), GeoCotrans.calculate_gst_alt(dt), rtol = 5.0e-4)
        # test csundir
        @test all(isapprox.(GeoCotrans.csundir_astrolib(dt), GeoCotrans.csundir(dt), rtol = 5.0e-4))
    end

    using Chairmarks
    dt = Date(2021, 3, 28)

    @info @b GeoCotrans.csundir_astrolib($dt)
    @info @b GeoCotrans.csundir($dt)
end

@testitem "IGRF get_B" begin
    using Dates

    # Test IGRF magnetic field calculation
    r, θ, φ = 6500.0, 30.0, 4.0
    t = Date(2021, 3, 28)
    B_true = (-46077.31133522, -14227.12618499, 233.14355744)
    @test all(igrf_Bd(r, θ, φ, t) .≈ B_true)

    𝐫 = GDZ(0, 60.39299, 5.32415)
    B_true = (458.89660058, 14996.72893889, -49019.55372591)
    @test all(igrf_B(𝐫, t) .≈ B_true)
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

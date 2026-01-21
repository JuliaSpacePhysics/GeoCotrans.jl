@testitem "IGRF coefficients" begin
    using Dates
    t = DateTime(2021, 3, 28, 1)
    g, h = get_igrf_coeffs(t)
    @test get_igrf_coeffs(t) != get_igrf_coeffs(Date(t))
    @test get_igrf_coeffs(DateTime(2021, 3, 28)) == get_igrf_coeffs(Date(t))
    @test g isa AbstractArray
    @test length(g) == 105
    @test length(h) == 105
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

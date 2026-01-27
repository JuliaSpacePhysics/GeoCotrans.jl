@testitem "get_dipole_terms" begin
    using Dates
    using GeoCotrans: get_dipole_terms

    # Test with known IGRF coefficients from a specific date
    t = DateTime(2015)
    g, h = get_igrf_coeffs(t)

    result = get_dipole_terms(g, h)

    # Verify b0 (dipole strength) is physically reasonable for Earth (~30,000 nT)
    @test 29800 < result.b0 < 30000
    # Verify offset terms (x0, y0, z0) are small compared to Earth's radius
    # These represent the dipole offset in Earth radii, should be << 1
    @test abs(result.x0) < 0.1
    @test abs(result.y0) < 0.1
    @test abs(result.z0) < 0.1

    # Test with a different epoch to ensure consistency
    t2 = Date(2010)
    g2, h2 = get_igrf_coeffs(t2)
    result2 = get_dipole_terms(g2, h2)
    # Dipole parameters should be similar across nearby epochs
    @test 29800 < result2.b0 < 30000

    # Cross validation: parameters should be consistent across nearby epochs
    @test result.θ ≈ result2.θ rtol = 1.0e-1
    @test result.φ ≈ result2.φ rtol = 1.0e-1
    @test result.x0 ≈ result2.x0 rtol = 1.0e-1
    @test result.y0 ≈ result2.y0 rtol = 1.0e-1
    @test result.z0 ≈ result2.z0 rtol = 1.0e-1
    @info result result2

end

@testitem "Coordinate transformations" begin
    using GeoCotrans: car2sph, car2sphd, sph2car
    # Test round-trip: Cartesian -> Spherical -> Cartesian
    𝐫 = [1.0, 2.0, 3.0]
    𝐫̂ = car2sph(𝐫)
    𝐫̂2 = sph2car(𝐫̂)
    @test 𝐫 ≈ 𝐫̂2 atol = 1.0e-10

    # Test specific cases
    # Point on z-axis
    (r, θ, φ) = car2sph(0.0, 0.0, 5.0)
    @test r ≈ 5.0
    @test θ ≈ 0.0  # Colatitude = 0 at north pole

    # Point on x-axis
    𝐫 = [3.0, 0.0, 0.0]
    (r, θ, φ) = car2sph(𝐫)
    @test r ≈ 3.0
    @test θ ≈ π / 2
    @test φ ≈ 0.0
    @test car2sphd(𝐫) == [3.0, 90, 0.0]

    # Point on y-axis
    (r, θ, φ) = car2sph(0.0, 4.0, 0.0)
    @test r ≈ 4.0
    @test θ ≈ π / 2
    @test φ ≈ π / 2

    # Origin
    (r, θ, φ) = car2sph(0.0, 0.0, 0.0)
    @test r ≈ 0.0
end


@testitem "Field transformations" begin
    using GeoCotrans: car2sph, car2sphd, sph2car

    # Test roundtrip at a known position
    θ, φ = π / 4, π / 6
    rθφ = [1.0, θ, φ]
    B = [100.0, 50.0, -30.0]

    # Spherical to Cartesian
    B_cart = sph2car(B, rθφ)

    # Get Cartesian position for reverse transform
    𝐱 = sph2car(1.0, θ, φ)
    # Cartesian to Spherical
    B_sph = car2sph(B_cart, 𝐱)

    @test B_sph ≈ B atol = 1.0e-10

    # polar axis: φ = 0, matching car2sph(x, y, z)
    for θ in (0.0, π)
        @test car2sph(sph2car(B, [2.0, θ, 0.0]), sph2car(2.0, θ, 0.0)) ≈ B
    end
end

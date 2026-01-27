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

@testitem "IGRF calculation" begin
    using Dates
    using GeoCotrans: IGRF, R🜨, igrf_B, igrf_Benu, igrf_Bd
    using GeoCotrans: Spherical, car2sph
    using Chairmarks

    # Test IGRF magnetic field calculation
    r = 1.0
    θ = deg2rad(45)
    φ = deg2rad(45)
    t = Date(2015)
    B = igrf_B(r, θ, φ, t)
    B_true = [-45469.44626375856, -21942.539310375545, 2933.49091800253]
    @test B ≈ B_true
    model = IGRF()
    @test model(r, θ, φ, t) ≈ B

    @info @b $igrf_B($r, $θ, $φ, $t), $model($r, $θ, $φ, $t)

    r, θ, φ = 1, 45, 45  # r in Re
    @test igrf_Bd(r, θ, φ, t) ≈ B_true

    𝐫 = GDZ(60.39299, 5.32415)
    t = Date(2021, 3, 28)
    B_true = [458.89660058, 14996.72893889, -49019.55372591]
    @test igrf_Benu(𝐫, t) ≈ B_true
    @test model(𝐫, t) ≈ B_true
    gdz_array = [𝐫;; 𝐫]
    @test model(gdz_array, [t, t]; in = GDZ()) ≈ [B_true B_true]
    # igrf_Benu(𝐫, t) bypass conversion using Cartesian so should be faster
    @info @b igrf_Benu($𝐫, $t), $model($𝐫, $t)

    # Test (GSM, Cartesian3) input
    t = DateTime("1970-01-01T00:01:40")
    r = GSM(1, 2, 3)
    B_true = [-548.8914669609419, -572.8115580345975, -319.8510173486891]
    @test GeoCotrans.igrf(r, t) ≈ B_true
    @test model(r, t) ≈ B_true
    @info @b $igrf($r, $t), $IGRF()($r, $t)

    # Test (GSM, Spherical) input
    rθφ = car2sph(r)
    in_csys = (GSM(), Spherical())
    @test model(rθφ, t; in = in_csys) == model(rθφ, t; in = in_csys, out = in_csys)
    @test model(rθφ, t; in = in_csys) ≈ car2sph(B_true, r)
end

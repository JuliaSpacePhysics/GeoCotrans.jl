@testitem "IGRF coefficients" begin
    using Dates
    t = DateTime(2021, 3, 28, 1)
    g, h = get_igrf_coeffs(t)
    @test get_igrf_coeffs(t) != get_igrf_coeffs(Date(t))
    @test get_igrf_coeffs(DateTime(2021, 3, 28)) == get_igrf_coeffs(Date(t))
    @test g isa AbstractArray
    @test length(g) == 105
    @test length(h) == 105
    g!, h! = zeros(105), zeros(105)
    get_igrf_coeffs!(g!, h!, t)
    @test (g!, h!) == (g, h)
end

@testitem "IGRF epochs" begin
    using Dates
    m = IGRF()
    # 2025–2030: IGRF-14 secular variation, g10 changes by 12.6 nT/yr
    @test get_igrf_coeffs(DateTime(2027, 7, 2, 12))[1][2] ≈ -29350.0 + 2.5 * 12.6 atol = 0.05
    @test m(1.0, 0.5, 0.5, DateTime(2030)) isa AbstractVector
    @test_throws DomainError m(1.0, 0.5, 0.5, DateTime(2030, 1, 2))
    @test_throws DomainError m(1.0, 0.5, 0.5, DateTime(1964, 12, 31))
end

@testitem "SphericalHarmonicModel" begin
    using Dates
    using GeoCotrans: SphericalHarmonicModel, GaussCoefficients, IGRF_COEFFS, igrf_ga, igrf_ha, igrf_B, get_coeff, R🜨
    g, h = igrf_ga[2015], igrf_ha[2015]
    G, H = ([m <= n ? get_coeff(c, n, m) : 0.0 for n in 0:13, m in 0:13] for c in (g, h))
    earth = (; radius = R🜨, frame = GEO())
    packed = SphericalHarmonicModel("IGRF2015", GaussCoefficients(g, h); earth...)
    matrix = SphericalHarmonicModel("IGRF2015", GaussCoefficients(G, H); earth...)
    @test packed.degree == packed.order == 13
    @test SphericalHarmonicModel("axisymmetric", GaussCoefficients(G[:, 1:1], H[:, 1:1]); earth...).order == 0
    @test packed.coeffs[4, 2] == matrix.coeffs[4, 2] == (g = G[5, 3], h = H[5, 3])
    @test_throws BoundsError matrix.coeffs[14, 0]
    @test_throws BoundsError packed.coeffs[1, 2]
    @test_throws BoundsError IGRF_COEFFS(DateTime(2020))[14, 0]
    @test_throws ArgumentError GaussCoefficients(g[2:end], h[2:end])
    B = igrf_B(1.0, 0.7, 2.0, Date(2015))
    @test packed(1.0, 0.7, 2.0) ≈ B
    @test matrix(1.0, 0.7, 2.0) ≈ B

    snapshot(y) = SphericalHarmonicModel("IGRF$y", GaussCoefficients(igrf_ga[y], igrf_ha[y]); earth...)(1.0, 0.7, 2.0)
    t0, t, t1 = DateTime(2025), DateTime(2027, 7, 2, 12), DateTime(2030)
    α = (t - t0) / (t1 - t0)
    @test IGRF()(1.0, 0.7, 2.0, t) ≈ (1 - α) * snapshot(2025) + α * snapshot(2030)
end

@testitem "Unitful positions" begin
    using Dates, Unitful
    using GeoCotrans: SphericalHarmonicModel, IGRF_COEFFS, R🜨
    t = DateTime(2020)
    r, θ, φ = 1.5, 0.7, 2.0
    x = r * R🜨 * u"km"
    m = IGRF()
    @test m(uconvert(u"m", x), θ, φ, t) ≈ m(r, θ, φ, t)
    @test m(uconvert(u"m", x), θ, φ, t; from = GEO()) ≈ m(r, θ, φ, t)
    @test m(x, 0u"km", 2x, t; from = GSM()) ≈ m(r, 0.0, 2r, t; from = GSM())
end

@testitem "evalsph is -∇V" begin
    using Dates
    using ForwardDiff
    using SatelliteToolboxLegendre: legendre
    using GeoCotrans: evalsph, get_coeff

    # Reference independent of evalsph's fused recursion: V from SatelliteToolboxLegendre, differentiated by AD
    function V(G, H, (r, θ, φ), n, mmax)
        P = legendre(Val(:schmidt), θ, n)
        return sum(
            (1 / r)^(l + 1) * P[l + 1, m + 1] * (get_coeff(G, l, m) * cos(m * φ) + get_coeff(H, l, m) * sin(m * φ))
                for l in 1:n for m in 0:min(l, mmax)
        )
    end
    function B_ref(G, H, r, θ, φ, n, mmax = n)
        dV = ForwardDiff.gradient(x -> V(G, H, x, n, mmax), [r, θ, φ])
        return [-dV[1], -dV[2] / r, -dV[3] / (r * sin(θ))]
    end

    g, h = get_igrf_coeffs(Date(2015))
    G, H = randn(31, 31), randn(31, 31)
    for r in (1.0, 2.5), θ in (1.0e-6, 0.3, π / 4, 2.0, π - 1.0e-6), φ in (0.0, 1.0, 4.0)
        @test evalsph((; g, h), r, θ, φ, 13) ≈ B_ref(g, h, r, θ, φ, 13)
        @test evalsph((; g = G, h = H), r, θ, φ, 30, 7) ≈ B_ref(G, H, r, θ, φ, 30, 7)
    end
    @test_throws BoundsError evalsph((; g, h), 1.0, 0.5, 0.5, 14)
    # At the poles B_θ, B_φ are the limits along the meridian φ
    for φ in (0.0, 1.0)
        @test evalsph((; g, h), 1.0, 0.0, φ, 13) ≈ evalsph((; g, h), 1.0, 1.0e-9, φ, 13)
        @test evalsph((; g, h), 1.0, π, φ, 13) ≈ evalsph((; g, h), 1.0, π - 1.0e-9, φ, 13)
    end
end

@testitem "IGRF calculation" begin
    using Dates
    using GeoCotrans: IGRF, R🜨, igrf_B, igrf_Benu, igrf_Bd
    using GeoCotrans: Spherical, car2sph, sph2car
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

    # Cartesian output
    @test model(r, θ, φ, t; to = :cartesian) ≈ sph2car(B, [r, θ, φ])
    @test (@b $model($r, $θ, $φ, $t; to = :cartesian)).allocs == 0

    r, θ, φ = 1, 45, 45  # r in Re
    @test igrf_Bd(r, θ, φ, t) ≈ B_true

    @testset "GDZ and Array input" begin
        𝐫 = GDZ(60.39299, 5.32415)
        t = Date(2021, 3, 28)
        B_true = [458.89660058, 14996.72893889, -49019.55372591]
        @test igrf_Benu(𝐫, t) ≈ B_true
        @test model(𝐫, t) ≈ B_true
        @test model(GDZ(90, 0, 100), t) ≈ igrf_Benu(GDZ(90, 0, 100), t)
        gdz_array = [𝐫;; 𝐫]
        @test model(gdz_array, [t, t]; from = GDZ()) ≈ [B_true B_true]
        @test model(gdz_array, t; from = GDZ()) ≈ [B_true B_true]
        # igrf_Benu(𝐫, t) bypass conversion using Cartesian so should be faster
        @info @b igrf_Benu($𝐫, $t), $model($𝐫, $t)
    end

    @testset "(GSM, Cartesian3) input" begin
        t = DateTime("1970-01-01T00:01:40")
        r = GSM(1, 2, 3)
        B_true = [-548.8914669609419, -572.8115580345975, -319.8510173486891]
        @test GeoCotrans.igrf(r, t) ≈ B_true
        @test model(r, t) ≈ B_true
        @info @b $igrf($r, $t), $IGRF()($r, $t)

        # Test (GSM, Spherical) input
        rθφ = car2sph(r)
        in_csys = (GSM(), Spherical())
        @test model(rθφ, t; from = in_csys) == model(rθφ, t; from = in_csys, to = in_csys)
        @test model(rθφ, t; from = in_csys) ≈ car2sph(B_true, r)
    end
end

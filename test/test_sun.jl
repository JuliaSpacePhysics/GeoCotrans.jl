@testitem "Sun direction and GMST vs AstroLib" begin
    using Dates
    using AstroLib
    using LinearAlgebra
    using GeoCotrans: sun_gei, gmst

    # Reference: https://aa.usno.navy.mil/faq/GAST, https://github.com/JuliaAstro/AstroLib.jl/blob/main/src/ct2lst.jl
    angle_deg(a, b) = rad2deg(atan(norm(a × b), a ⋅ b))
    for t in DateTime(1965):Hour(1001):DateTime(2100)
        jd = jdcnv(t)
        ra, dec = sunpos(jd; radians = true)
        @test angle_deg(sun_gei(t), [cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)]) < 0.015
        @test rad2deg(abs(rem2pi(gmst(t) - ct2lst(0, jd) * 2π / 24, RoundNearest))) < 0.003
    end
end

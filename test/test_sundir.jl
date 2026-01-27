@testitem "Validation with AstroLib" begin
    using Dates
    using AstroLib

    # Reference: https://aa.usno.navy.mil/faq/GAST, https://github.com/JuliaAstro/AstroLib.jl/blob/main/src/ct2lst.jl
    function calculate_gst(time)
        # ct2lst returns local sidereal time in hours (0-24)
        # For Greenwich, longitude=0, so GST = LST
        return ct2lst(0, jdcnv(time)) * 2π / 24
    end

    function csundir_astrolib(time)
        jd = jdcnv(time)
        gst = ct2lst(0, jd) * 2π / 24
        return gst, sunpos(jd; radians = true)...
    end

    for yr in 2010:2025
        dt = Date(yr, 1, 1)
        # test GST
        @test isapprox(calculate_gst(dt), GeoCotrans.calculate_gst_alt(dt), rtol = 5.0e-4)
        # test csundir
        @test all(isapprox.(csundir_astrolib(dt), GeoCotrans.csundir(dt), rtol = 5.0e-4))
    end

    using Chairmarks
    dt = Date(2021, 3, 28)

    @info @b GeoCotrans.csundir($dt), csundir_astrolib($dt)
end

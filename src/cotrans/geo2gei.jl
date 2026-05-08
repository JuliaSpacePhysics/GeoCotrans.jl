# https://pyspedas.readthedocs.io/en/latest/coords.html#pyspedas.cotrans_tools.cotrans_lib.subgeo2gei

"""
    rotation(GEI, GEO, gst_or_time)

GEI → GEO rotation. Accepts either Greenwich sidereal time `gst` (radians) or
an `AbstractTime` from which `gst` is derived.
"""
@inline function rotation(::Type{GEI}, ::Type{GEO}, gst)
    sgst, cgst = sincos(gst)
    return SA[cgst sgst 0; -sgst cgst 0; 0 0 1]
end

rotation(::Type{GEI}, ::Type{GEO}, time::AbstractTime) = rotation(GEI, GEO, calculate_gst_alt(time))

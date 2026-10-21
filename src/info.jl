const coord_text = Dict(
    :geo => "Geographic (GEO)",
    :gei => "Geocentric Equatorial Inertial (GEI)",
    :gse => "Geocentric Solar Ecliptic (GSE)",
    :gsm => "Geocentric Solar Magnetospheric (GSM)",
    :mag => "Geomagnetic (MAG)",
    :sm => "Solar Magnetic (SM)",
    :j2000 => "J2000 Geocentric Equatorial (J2000)",
)

const FrameDescriptions = Dict(
    :GEI => "Geocentric Equatorial Inertial (GEI) reference frame.",
    :GEO => "Geocentric Geographic (GEO) reference frame.",
    :GSM => """Geocentric Solar Magnetospheric (GSM) reference frame.""",
    :GSE => "Geocentric Solar Ecliptic (GSE) reference frame.",
    :MAG => "Geomagnetic (MAG) reference frame.",
    :SM => "Solar Magnetic (SM) reference frame.",
    :J2000 => "J2000 Geocentric Equatorial (J2000, EME2000) reference frame.",
)

const FrameDefinitions = Dict(
    :GEI => "Z-axis is parallel to Earth's rotation axis (positive northward). X-axis points toward the vernal equinox (first point of Aries) of date. As in SPEDAS, mean and true of date (differing by nutation, < 20″, under the 0.01° accuracy) are not distinguished; J2000 maps to true of date.",
    :J2000 => "Like GEI, but fixed at the mean equator and equinox of the J2000.0 epoch (2000-01-01T12:00 TT). Inertial.",
    :GEO => "Z-axis is parallel to Earth's rotation axis (positive northward). X-axis passes through the intersection of the equator and the Greenwich meridian. Rotates with Earth.",
    :GSE => "X points sunward from Earth's center. Z-axis is perpendicular to the ecliptic plane (positive north).",
    :SM => "Z-axis is parallel to Earth's magnetic dipole axis (positive northward). The X-Z plane contains the Earth-Sun direction.",
    :GSM => "X points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).",
    :MAG => "Z-axis is parallel to Earth's magnetic dipole axis (positive northward). Y-axis is perpendicular to the plane containing the dipole and Earth's rotation axis."
)

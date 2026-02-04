const coord_text = Dict(
    :geo => "Geographic (GEO)",
    :gei => "Geocentric Equatorial Inertial (GEI)",
    :gse => "Geocentric Solar Ecliptic (GSE)",
    :gsm => "Geocentric Solar Magnetospheric (GSM)",
    :mag => "Geomagnetic (MAG)",
    :sm => "Solar Magnetic (SM)"
)

const FrameDescriptions = Dict(
    :GEI => "Geocentric Equatorial Inertial (GEI) reference frame.",
    :GEO => "Geocentric Geographic (GEO) reference frame.",
    :GSM => """Geocentric Solar Magnetospheric (GSM) reference frame.""",
    :GSE => "Geocentric Solar Ecliptic (GSE) reference frame.",
    :MAG => "Geomagnetic (MAG) reference frame.",
    :SM => "Solar Magnetic (SM) reference frame.",
)

const FrameDefinitions = Dict(
    :GSM => "X points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).",
    :MAG => "Z-axis is parallel to Earth's magnetic dipole axis (positive northward). Y-axis is perpendicular to the plane containing the dipole and Earth's rotation axis."
)

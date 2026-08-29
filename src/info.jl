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
    :GEI => "Z-axis is parallel to Earth's rotation axis (positive northward). X-axis points toward the vernal equinox (first point of Aries). Non-rotating.",
    :GEO => "Z-axis is parallel to Earth's rotation axis (positive northward). X-axis passes through the intersection of the equator and the Greenwich meridian. Rotates with Earth.",
    :GSE => "X points sunward from Earth's center. Z-axis is perpendicular to the ecliptic plane (positive north).",
    :SM => "Z-axis is parallel to Earth's magnetic dipole axis (positive northward). The X-Z plane contains the Earth-Sun direction.",
    :GSM => "X points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).",
    :MAG => "Z-axis is parallel to Earth's magnetic dipole axis (positive northward). Y-axis is perpendicular to the plane containing the dipole and Earth's rotation axis."
)

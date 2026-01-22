"""
    igrf(𝐫, t)

Calculate IGRF model given coordinates `𝐫` at time `t`.
Output magnetic field in the frame of `𝐫`.

    igrf(𝐫::CoordinateVector{GSM}, t) -> (Bx, By, Bz) in GSM frame
    igrf(𝐫::CoordinateVector{SPH}, t) -> (Br, Bθ, Bφ) in GEO frame
    igrf(𝐫::CoordinateVector{GDZ}, t) -> (Be, Bn, Bu) in ENU frame

GDZ: Geodetic
ENU: East-North-Up
"""
function igrf end
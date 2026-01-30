using LinearAlgebra
using Bumper
using SatelliteToolboxLegendre: legendre!

"""
    evalsph(coeffs, r, θ, φ)

Evaluate the magnetic field at a point in spherical coordinates.

# Arguments
- `coeffs`: The magnetic field coefficients
- `r`: Radial distance in planetary radii (dimensionless)
- `θ`: Colatitude in radians [0, π]
- `φ`: East longitude in radians [0, 2π]

# Returns
- `SVector{3,Float64}`: Magnetic field vector [B_r, B_θ, B_φ] in nanoTesla

# Coordinate System
- Spherical coordinates (r, θ, φ) where:
  - r: radial distance from planet center
  - θ: colatitude (0 at north pole, π at south pole)
  - φ: east longitude (0 at prime meridian, increases eastward)

# Mathematical Formulation
The magnetic field is derived from the scalar potential:
```math
B_r = -\\frac{∂V}{∂r}, \\quad
B_θ = -\\frac{1}{r}\\frac{∂V}{∂θ}, \\quad
B_φ = -\\frac{1}{r\\sin θ}\\frac{∂V}{∂φ}
```
"""
@inline function evalsph(coeffs, r, θ, φ, max_degree, max_order = max_degree)
    T = promote_type(eltype(r), eltype(θ), eltype(φ))
    return @no_escape begin
        P = @alloc(T, max_degree + 1, max_degree + 1)
        sincos_mφs = @alloc(Tuple{T, T}, max_degree + 1)
        evalsph!(P, sincos_mφs, coeffs.g, coeffs.h, promote(r, θ, φ)..., max_degree, max_order)
    end
end


# Accessor for 2D arrays
@inline get_coeff(G::AbstractMatrix, l, m) = G[l + 1, m + 1]

# Accessor for 1D vectors with triangular indexing: k = l*(l+1)÷2 + 1 + m
@inline get_coeff(G::AbstractVector, l, m) = G[l * (l + 1) ÷ 2 + 1 + m]

@inline function evalsph!(P, sincos_mφs, G, H, r::T, θ::T, φ::T, max_degree, max_order) where {T}
    sinθ, cosθ = sincos(θ)
    for m in eachindex(sincos_mφs)
        sincos_mφs[m] = sincos((m - 1) * φ)
    end
    legendre!(Val(:schmidt), P, θ, max_degree)

    Br, Bθ, Bφ = 0.0, 0.0, 0.0

    # r is in planetary radii (dimensionless), so (1/r)^k gives the radial dependence
    # For internal field: B components use (1/r)^(n+2)
    ratio = 1.0 / r
    pow = ratio * ratio * ratio   # (1/r)^(n+2) for magnetic field
    flag = abs(sinθ) > 1.0e-10
    sinθ = flag ? sinθ : 1.0e-10
    # Sum over degrees and orders
    @inbounds for l in 1:max_degree
        Vl, dVl, φVl = 0.0, 0.0, 0.0
        for m in 0:min(l, max_order)
            g = get_coeff(G, l, m)
            h = get_coeff(H, l, m)
            sin_mφ, cos_mφ = sincos_mφs[m + 1]
            Pₗₘ = P[l + 1, m + 1]

            # Spherical harmonic term
            Y = g * cos_mφ + h * sin_mφ
            # Derivatives
            dY_dφ = m * (-g * sin_mφ + h * cos_mφ)

            dPₗₘ = if m == l
                l * cosθ * Pₗₘ / sinθ
            else
                Pₗ₋₁ₘ = P[l, m + 1]
                s = sqrt((l + m) * (l - m))
                (l * cosθ * Pₗₘ - s * Pₗ₋₁ₘ) / sinθ
            end

            # Contribution to field components
            # B_r = -∂V/∂r: factor of (n+1) from derivative of (a/r)^(n+1)
            Vl += Y * Pₗₘ
            # B_θ = -(1/r)∂V/∂θ
            dVl += Y * dPₗₘ
            # B_φ = -(1/(r sin θ))∂V/∂φ
            flag && (φVl -= dY_dφ * Pₗₘ)
        end
        Br += (l + 1) * pow * Vl
        Bθ -= pow * dVl
        Bφ += pow * φVl
        pow *= ratio
    end
    return SVector{3, T}(Br, Bθ, flag ? Bφ / sinθ : Bφ)
end

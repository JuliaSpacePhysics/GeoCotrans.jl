# GeoCotrans.jl Use Case Workflow

Guide for common use cases: coordinate transformations, IGRF field evaluation, and field line tracing.

## Setup

```julia
using GeoCotrans
using Dates
```

## Coordinate Transformations

Convert a position vector between geocentric frames using a `DateTime` for time-dependent rotations:

```julia
t = DateTime(2020, 6, 1)

# GEO → GSM
pos_geo = CoordinateVector{GEO}(3.0, 1.0, 2.0, t)   # RE units
pos_gsm = geo2gsm(pos_geo)

# GSM → SM
pos_sm = gsm2sm(pos_gsm)

# Chain: GEO → GEI → GSE
pos_gse = geo2gse(pos_geo)
```

Available frames: `GEO`, `GEI`, `GSM`, `GSE`, `MAG`, `SM`

Available representations: `Cartesian3` (default), `Spherical`, `Geodetic`

```julia
# Spherical representation
pos_sph = CoordinateVector{GEO, Spherical}(6371.0, 45.0, 90.0, t)  # r[km], θ[deg], φ[deg]
```

## IGRF Field Evaluation

Evaluate the IGRF-14 geomagnetic field at a given location and time (valid 1965–2030):

```julia
# igrf(date, altitude_km, colat_deg, lon_deg) → (Br, Bθ, Bφ) in nT
Br, Bθ, Bφ = igrf(DateTime(2020, 1, 1), 400.0, 45.0, 120.0)

# Or use the IGRF field model object
field = IGRF()
B = field(pos_geo)   # returns field vector in same frame
```

## Field Line Tracing

Requires SciML (OrdinaryDiffEq) loaded as an extension:

```julia
using OrdinaryDiffEq

# Trace a field line from a starting position
start = CoordinateVector{GEO}(2.0, 0.0, 0.0, t)   # 2 RE on equator
trace = trace_fieldline(start, IGRF(); alg=Tsit5())

# Result is a trajectory of CoordinateVectors
```

## DimensionalData Integration

Apply transformations to labeled arrays via the DimensionalData extension:

```julia
using DimensionalData

# Transform a DimArray of positions
da = DimArray(positions, (Ti(times),))
da_gsm = geo2gsm(da)
```

## Common Patterns

- All transformation functions accept both `CoordinateVector` and plain `SVector`/array inputs
- Time is embedded in `CoordinateVector` and propagated automatically through chains
- Inverse transforms are available: `gsm2geo`, `gsm2gei`, etc.
- Use `mlt(pos)` to compute Magnetic Local Time from a position

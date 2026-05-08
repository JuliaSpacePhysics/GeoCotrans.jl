# GeoCotrans.jl

Julia package for transformations between geocentric coordinate systems (GEI/GEO/GSE/GSM/MAG/SM), IGRF-14 geomagnetic field evaluation, and magnetic field line tracing.

## Quick Start

```julia
using Pkg; Pkg.add("GeoCotrans")
using GeoCotrans, Dates

times = Date.(2015:2020)
x_geo = rand(3, 6)

x_gsm = geo2gsm(x_geo, times) # equivalent to transform(GSM, GEO, x_geo, times)
get_mlt(x_geo, times) # Magnetic Local Time

B_gsm = igrf(x_gsm, times; in = GSM())

using OrdinaryDiffEqTsit5, CairoMakie

sol = trace(GEO(3.0, 0.0, 0.0), DateTime(2020, 1, 1), Tsit5()) # Trace from [3, 0, 0] Earth radii
plot(sol; idxs = (1, 2)) # Equatorial plane (X-Y)
```

## API Map

- Reference Frames Types:
  - `GEI`, `GEO`, `GSE`, `GSM`, `MAG`, `SM`, `GDZ`
- Coordinate transforms `transform(to, from, pos, time)`
  - Utility functions (each has inverse `b2a`)
    - `geo2gei`, `geo2gsm`, `geo2mag`
    - `gei2gsm`, `gse2gsm`, `gsm2sm`
    - `gei2sm`, `geo2sm`, `gei2mag`

See [documentation](https://JuliaSpacePhysics.github.io/GeoCotrans.jl/dev/) for full signatures.

## Conventions

- GDZ(lat, lon, alt) uses geodetic latitude and longitude in degrees, altitude in kilometers.
- Other inputs use radius in Earth radii and angles in radians.

## Reference

- [SSC Appendix C — Description of Coordinate Systems](https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml)
- [Coordinate transformations between geocentric systems](https://www.mssl.ucl.ac.uk/grid/iau/extra/local_copy/SP_coords/geo_tran.htm)
- [SPEDAS.jl — Coordinates explanation](https://beforerr.github.io/SPEDAS.jl/dev/explanations/coords/)
- [IAGA / NOAA NCEI — IGRF](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field)

## Related packages

- [SatelliteToolboxGeomagneticField.jl](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl) — IGRF-13 and dipole models (Julia)
- [ppigrf](https://github.com/IAGA-VMOD/ppigrf) — pure Python IGRF
- [geopack](https://github.com/tsssss/geopack) — Python IGRF + Tsyganenko
- [TREPS](https://cdpp.irap.omp.eu/index.php/services/treps) — online SPICE-based transforms
- [pymaginverse](https://github.com/outfrenk/pymaginverse) — geomagnetic field inversion

## Status

[![DOI](https://zenodo.org/badge/1004669370.svg)](https://doi.org/10.5281/zenodo.15709873)
[![version](https://juliahub.com/docs/General/GeoCotrans/stable/version.svg)](https://juliahub.com/ui/Packages/General/GeoCotrans)
[![Build Status](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/GeoCotrans.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/GeoCotrans.jl)

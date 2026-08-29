# GeoCotrans.jl

[![DOI](https://zenodo.org/badge/1004669370.svg)](https://doi.org/10.5281/zenodo.15709873)

Julia package for transformations between geocentric coordinate systems (GEI/GEO/GSE/GSM/MAG/SM), IGRF-14 geomagnetic field evaluation, and magnetic field line tracing.

## Quick Start

```julia
using Pkg; Pkg.add("GeoCotrans")
using GeoCotrans, Dates

times = Date.(2015:2020)
x_geo = rand(3, 6)

x_gsm = geo2gsm(x_geo, times) # equivalent to transform(GEO() => GSM(), x_geo, times)
get_mlt(x_geo, times) # Magnetic Local Time

model = IGRF()
B_gsm = model(x_gsm, times; from = GSM()) # nT

using TsyganenkoModels
model = TsyIGRF(T89(2)) # Kp level = 2

using OrdinaryDiffEqTsit5, CairoMakie

sol = trace(GEO(3.0, 0.0, 0.0), DateTime(2020, 1, 1), Tsit5()) # Trace from [3, 0, 0] Earth radii
plot(sol; idxs = (1, 2)) # Equatorial plane (X-Y)
```

## API

Full API map with units and conventions: `@doc(GeoCotrans)`, or the [documentation](https://JuliaSpacePhysics.github.io/GeoCotrans.jl/dev/) (including [validation](https://juliaspacephysics.github.io/GeoCotrans.jl/dev/coords/) against IRBEM and PySPEDAS).

## Reference

- [SSC Appendix C — Description of Coordinate Systems](https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml)
- [SPEDAS.jl — Coordinates explanation](https://beforerr.github.io/SPEDAS.jl/dev/explanations/coords/)

## Related packages

- [TsyganenkoModels.jl](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl) for external Tsyganenko magnetosphere models
- [GeoAACGM.jl](https://github.com/JuliaSpacePhysics/GeoAACGM.jl) for AACGM coordinates
- [SatelliteToolboxGeomagneticField.jl](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl) — IGRF-13 and dipole models (Julia)
- [ppigrf](https://github.com/IAGA-VMOD/ppigrf) — pure Python IGRF
- [geopack](https://github.com/tsssss/geopack) — Python IGRF + Tsyganenko
- [TREPS](https://cdpp.irap.omp.eu/index.php/services/treps) — online SPICE-based transforms
- [pymaginverse](https://github.com/outfrenk/pymaginverse) — geomagnetic field inversion

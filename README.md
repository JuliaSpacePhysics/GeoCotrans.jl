# GeoCotrans.jl

[![DOI](https://zenodo.org/badge/1004669370.svg)](https://doi.org/10.5281/zenodo.15709873)
[![version](https://juliahub.com/docs/General/GeoCotrans/stable/version.svg)](https://juliahub.com/ui/Packages/General/GeoCotrans)

[![Build Status](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/GeoCotrans.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/GeoCotrans.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

A Julia package for coordinate transformations between common geocentric systems. It also provides interpolation of [International Geomagnetic Reference Field (IGRF)](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field) field model.

For more information, see the [Documentation](https://juliaspacephysics.github.io/GeoCotrans.jl/dev/) and
[Coordinates - SPEDAS](https://beforerr.github.io/SPEDAS.jl/dev/explanations/coords/).

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("GeoCotrans")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://JuliaSpacePhysics.github.io/GeoCotrans.jl/dev/)

## Reference

- [SSC: APPENDIX C: Description of Selected Coordinate Systems Used in SSC Programs](https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.shtml)

## Elsewhere

- [SatelliteToolboxGeomagneticField.jl](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl): Models to compute the geomagnetic field (IGRF-13, dipole model)
- [ppigrf](https://github.com/IAGA-VMOD/ppigrf): Pure Python code to calculate IGRF model predictions.
- [geopack](https://github.com/tsssss/geopack): Python code to calculate IGRF model predictions.
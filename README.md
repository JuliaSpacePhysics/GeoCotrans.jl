# GeoCotrans.jl

[![Build Status](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/GeoCotrans.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaspacephysics.github.io/GeoCotrans.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/GeoCotrans.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/GeoCotrans.jl)

A Julia package for coordinate transformations between common geocentric systems. It also provides interpolation of [International Geomagnetic Reference Field (IGRF)](https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field) field model.

For more information, see the [Documentation](https://juliaspacephysics.github.io/GeoCotrans.jl/dev/) and
[Coordinates - SPEDAS](https://beforerr.github.io/SPEDAS.jl/dev/explanations/coords/).

## Installation

```julia
using Pkg
Pkg.add("GeoCotrans")
```

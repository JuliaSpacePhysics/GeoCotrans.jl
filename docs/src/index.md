# GeoCotrans.jl

[![DOI](https://zenodo.org/badge/1004669370.svg)](https://doi.org/10.5281/zenodo.15709873)
[![version](https://juliahub.com/docs/General/GeoCotrans/stable/version.svg)](https://juliahub.com/ui/Packages/General/GeoCotrans)

```@docs
GeoCotrans
```

## Installation

```julia
using Pkg
Pkg.add("GeoCotrans")
```

## Usage

Magnetic Local Time (MLT) calculation

```@repl
using GeoCotrans, Dates

xGEO = [1., 2., 3.];
time = Date(2020);
get_mlt(xGEO, time)
```

## Coordinate Systems

```@docs; canonical=false
GEI
GEO
GSM
GSE
MAG
SM
```

```@docs; canonical=false
GDZ
```

## API Reference

### Public

```@autodocs
Modules = [GeoCotrans, GeoCotrans.FieldLineTracing]
Private = false
```

### Private

```@autodocs
Modules = [GeoCotrans]
Public = false
Order = [:function]
```

using BenchmarkTools
using Dates
using DimensionalData
using GeoCotrans

const SUITE = BenchmarkGroup()

const T0 = DateTime(2021, 3, 28, 1)
const TIMES = DateTime(2021, 1, 1, 0, 0, 0) .+ Hour.(0:127)

const GEO_POS = GEO(1.0, 2.0, 3.0, T0)
const GSM_POS = GSM(1.0, 2.0, 3.0, T0)
const R_SPH = [1.0, deg2rad(45), deg2rad(45)]
const R_GDZ = GDZ(60.39299, 5.32415)
const POS_TIME_MAJOR = rand(length(TIMES), 3)
const POS_COMPONENT_MAJOR = permutedims(POS_TIME_MAJOR)
const POS_DIMARRAY = DimArray(POS_TIME_MAJOR, (Ti(TIMES), Y(1:3)))

SUITE["transforms"] = BenchmarkGroup()
SUITE["transforms"]["scalar-direct"] = @benchmarkable gsm2sm($GSM_POS, $T0)
SUITE["transforms"]["scalar-chain"] = @benchmarkable geo2gsm($GEO_POS, $T0)
SUITE["transforms"]["batch-time-major"] = @benchmarkable geo2gsm($POS_TIME_MAJOR, $TIMES; dims = 1)
SUITE["transforms"]["batch-component-major"] = @benchmarkable geo2gsm($POS_COMPONENT_MAJOR, $TIMES; dims = 2)
SUITE["transforms"]["dimensionaldata"] = @benchmarkable geo2gsm($POS_DIMARRAY)

SUITE["igrf"] = BenchmarkGroup()
SUITE["igrf"]["scalar-spherical"] = @benchmarkable igrf($R_SPH, $T0)
SUITE["igrf"]["scalar-gdz"] = @benchmarkable igrf($R_GDZ, $T0)
SUITE["igrf"]["batch-time-major"] = @benchmarkable igrf($POS_TIME_MAJOR, $TIMES; in = (GEO(), Cartesian3()), dim = 1)
SUITE["igrf"]["batch-component-major"] = @benchmarkable igrf($POS_COMPONENT_MAJOR, $TIMES; in = (GEO(), Cartesian3()), dim = 2)

SUITE["mlt"] = BenchmarkGroup()
SUITE["mlt"]["scalar"] = @benchmarkable get_mlt($GEO_POS, $T0)
SUITE["mlt"]["batch-time-major"] = @benchmarkable get_mlt($POS_TIME_MAJOR, $TIMES; dim = 1)
SUITE["mlt"]["dimensionaldata"] = @benchmarkable get_mlt($POS_DIMARRAY)
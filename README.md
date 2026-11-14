# TsunamiTools.jl

[![CI](https://github.com/CyprienBosserelle/Tsunami.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/CyprienBosserelle/Tsunami.jl/actions/workflows/CI.yml)

Julia package with tools for tsunami science, including generating the initial
tsunami wave from an earthquake fault with the Okada (1985) solution.

## Installation

```julia
using Pkg
Pkg.add("TsunamiTools")
```

## Quick demo

```julia
using TsunamiTools
tsunamidemo()   # writes Vertical_displacement*.nc in the current directory
```

## Generating a tsunami initial wave

```julia
using TsunamiTools
using NetCDF

# Read bathy grid for a nearfield simulation
x  = NetCDF.ncread("Tohoku_bathy.nc", "lon")
y  = NetCDF.ncread("Tohoku_bathy.nc", "lat")
zb = NetCDF.ncread("Tohoku_bathy.nc", "z")

# Convert topography to water depth assuming 0.0 mean sea level
zs = 0.0
H = max.(zs, zs .- zb)

# Fault parameters relative to the bottom-left reference point
# faultparam(lon,lat,length,width,depth,strike,dip,rake,slip,tinit,trise)
# Note: bottom left corner of a fault with a strike of 192 means north-east corner.
# Slip is set to 0.0 here and computed from the magnitude below.
fault = faultparam(144.33, 39.6, 300.0, 150.0, 0.0, 192.0, 12.0, 90.0, 0.0, 0.0, 0.0)

# Convert length, width and depth from km to m
faultkm2m!(fault)

# Move from the bottom-left reference point to the centroid (what okada85 expects)
mvBLref2centroid!(fault)

# Slip for a given moment magnitude
Calcslip!(fault, 9.1)

# Initial tsunami wave (geographic coordinates; use InitTsunami for projected ones)
dz = InitTsunamiGeo(x, y, H, fault)

write2nc(x, y, dz, "Tohoku_fulldz.nc")
```

`Tohoku_fulldz.nc` can now be used to initialise a tsunami wave in a hydrodynamic model.

## Running the tests

```julia
using Pkg
Pkg.test("TsunamiTools")
```

## License

BSD 2-Clause, see [LICENSE](LICENSE). `src/Okada.jl` is a Julia translation of
François Beauducel's `okada85.m` (BSD License).

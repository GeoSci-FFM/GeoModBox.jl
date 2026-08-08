# Momentum Conservation Equation

This directory contains examples for solving the *Stokes equation* in one and two dimensions. The examples focus on pure Stokes flow and include, beyond others:

- [One-dimensional channel flow with constant and depth-dependent viscosity](./1D/ChannelFlow_1D.jl)
- [Two-dimensional falling block benchmark](./2D/FallingBlockBenchmark.jl)
- [Two-dimensional falling block example with constant viscosity using the defect correction method](./2D/FallingBlockConstEta_DC.jl)
- [Two-dimensional falling block example with variable viscosity using the defect correction method](./2D/FallingBlockVarEta_DC.jl)
- [Two-dimensional viscous inclusion problem](./2D/ViscousInclusion.jl)
- [Two-dimensional Rayleigh–Taylor instability benchmark](./2D/RTI.jl)
- [Van Keken Benchmark](./2D/VanKekenBenchmark.jl)

For detailed information on solving the *Stokes equation* using different numerical methods, please refer to the [GeoModBox.jl documentation](https://geosci-ffm.github.io/GeoModBox.jl/).
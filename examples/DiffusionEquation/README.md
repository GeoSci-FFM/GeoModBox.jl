# Energy Conservation Equation

This directory contains examples for solving the *diffusive component* of the *temperature conservation equation* in one and two dimensions. The included examples are, beyond others:

- Determination of one-dimensional geotherm profiles for [oceanic](./1D/OceanicGeotherm_1D.jl) and [continental](./1D/ContinentalGeotherm_1D.jl) settings.
- [Comparison of different finite difference (**FD**) schemes applied to a one-dimensional Gaussian temperature anomaly](./1D/Heat_1D_discretization.jl).
- [Backward](./2D/BackwardEuler.jl) and [forward](./2D/ForwardEuler.jl) Euler solutions using the defect correction method.
- [Two-dimensional resolution test for each **FD** scheme using a Gaussian temperature anomaly](./2D/Gaussian_Diffusion.jl).
- [Solution of a two-dimensional Poisson problem assuming variable thermal conductivity](./2D/Poisson_variable_k.jl).
- [Resolution test for the two-dimensional Poisson problem](./2D/Poisson_ResTest.jl).

For further details on the discretization methods and implementation of each example, please refer to the [GeoModBox.jl documentation](https://geosci-ffm.github.io/GeoModBox.jl/).

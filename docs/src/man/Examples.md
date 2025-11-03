# Examples and Benchmarks

The `GeoModBox.jl` provides various one- and two-dimensional examples and benchmark problems for each of the governing equations. The examples demonstrate how to implement different numerical solvers, apply scaling, and evaluate the advantages and limitations of various finite difference schemes.

By clicking on the title of each document page, you will be directed to the corresponding Julia script in the [examples directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples).

**Advection** 
- [2-D advection with constant velocity field](./examples/Advection2D.md)
- [Resolution test of 2-D advection](./examples/AdvectionRestest2D.md)

**Heat Diffusion** 
- [1-D continental geotherm](./examples/ContinentalGeotherm.md)
- [Comparison of FD schemes on a Gaussian anomaly](./examples/GaussianDiffusion1D.md)
- [1-D oceanic geotherm](./examples/OceanicGeotherm.md)
- [2-D Gaussian anomaly using iterative backward Euler solver](./examples/BackwardEuler_DC.md)
- [2-D Gaussian anomaly using iterative forward Euler solver](./examples/ForwardEuler_DC.md)
- [2-D resolution test with Gaussian anomaly](./examples/GaussianDiffusion2D.md)
- [2-D Poisson equation resolution test](./examples/PoissonRestest.md)
- [](./examples/PoissonVariablek.md)

**Thermal Convection Models** 
- [Bottom heated, isoviscous thermal convection](./examples/BottomHeatedConvection.md)
- [Internally heated, isoviscous thermal convection](./examples/InternallyHeatedConvection.md)
- [Mixed heated, isoviscous thermal convection](./examples/MixedHeatedConvection.md)

**Stokes Equation**
- [1D channel flow with constant and depth-dependent viscosity](./examples/ChannelFlow1D.md)
- [2D falling block benchmark](./examples/FallingBlockBenchmark.md)
- [2D falling block with constant or variable viscosity (defect correction)](./examples/FallingBlockDC.md)
- [2D Rayleigh–Taylor instability](./examples/RTI.md)
- [2D Rayleigh–Taylor instability benchmark](./examples/RTI_growth_rate.md)
- [2D viscous inclusion problem](./examples/ViscousInclusion.md)

In the following, the runtime for each of the provided examples is listed as a reference. 

| Example                           | Total Runtime                                         |
| :-------------------------------- | :---------------------------------------------------- |
| **Advection ===**                                                                         |
| 2D_Advection.jl                   | 1) Upwind: 5.88 s                                     |
|                                   | 2) SLF: 6.11 s                                        |
|                                   | 3) Semi-lag: 15.7 s                                   |
|                                   | 4) Tracers: 172 s                                     |
| 2D_Advection_ResolutionTest.jl    | 335 s                                                 |
| **Heat Diffusion ===**                                                                    |
| *1D* ---                                                                                  |
| ContinentalGeotherm_1D.jl         | 530 ms                                                | 
| Heat_1D_discretization.jl         | 2.28 s                                                |
| OceanicGeotherm_1D.jl             | 431 ms                                                |
| *2D* ---                                                                                  |
| BackwardEuler.jl                  | 29.9 s                                                | 
| ForwardEuler.jl                   | 5.89 s                                                |
| Gaussian_Diffusion.jl             | 678 s                                                 | 
| Poisson_RestTest.jl               | 34.7 s                                                |
| Poisson_variable_k.jl             | 5.24 s                                                |
| **Thermal Convection Models ===**                                                         |
| BottomHeated.jl                   | 6.37 h                                                |
| InternallyHeated.jl               | 5.94 h                                                |
| MixedHeated.jl                    |                                                       |
| **Stokes Equation ===**                                                                   |
| *1D* ---                                                                                  |
| ChannelFlow_1D.jl                 | 3.61 s                                                |
| *2D* ---                                                                                  |
| FallingBlockBenchmark.jl          | 1) Steady State: 8 s                                  |
|                                   | 2) Time-Dependent                                     |
|                                   |     a) Upwind: 74.4 s                                 |
|                                   |     b) SLF: 170 s                                     |
|                                   |     c) Semi-lag: 175 s                                |
|                                   |     d) Tracers: 233 s                                 |
| FallingBlockConstEta_DC.jl        | 332 ms                                                | 
| FallingBlockVarEta_DC.jl          | 25.8 s                                                |
| RTI.jl                            | 154 s                                                 |
| RTI_GrowthRate.jl                 | 49.4 s                                                |
| RTI_Growth_Rate_Res_Test.jl       | 250 s                                                 |
| RTI_Growth_Rate_Res_Test_2.jl     | 261 s                                                 |
| RTI_Growth_Rate_Res_Test_3.jl     | 446 s                                                 |

> **Note:** In `GeoModBox.jl`, thermal and kinematic boundary conditions are explicitly implemented within the solvers. The absolute values at *ghost nodes* are computed based on the values provided in the boundary condition tuple `BC`. Each tuple specifies the `type` (either Dirichlet or Neumann) and the corresponding `val`ue at each boundary.  
> For constant velocity boundary conditions, additional values must be defined in `val` for the boundary nodes (e.g., `BC.val.vxW`, `BC.val.vxE`). These additional values are required to **directly** solve the momentum equation using a the backslash operator and to update the right-hand side of the linear system. Furthermore, if non-zero, these values must be assigned to the initial boundary nodes of the respective velocity fields.  
> For more details on the implementation of constant velocity boundaries, refer to the documentation of the [viscous inclusion](./examples/ViscousInclusion.md) example or the [initial velocity setup](Ini.md).

> **Note:** By default, the results of time-dependent examples in `GeoModBox.jl` are stored as GIF animations. To visualize solutions at specific time steps without generating a GIF, set the parameter `save_fig = 0`. In this case, individual plots are not saved, so caution is advised when running problems that require multiple time step iterations.

> **Note:** Some examples use *named tuples* to define constants and parameters. Alternatively, *mutable structures* can be used—particularly useful when parameters need to be modified after initialization (e.g., for scaling purposes). A full transition from *named tuples* to *mutable structures* is planned for future versions of `GeoModBox.jl`.
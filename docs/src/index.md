# GeoModBox.jl

The **Geod**ynamic **Mod**elling Tool**Box** is a Julia package primarily intended for teaching purposes. It provides various finite difference, staggered discretization schemes to numerically solve the governing equations of two-dimensional geodynamic problems. These include the conservation equations of:

1) [**Energy**](./man/DiffMain.md), 
2) [**Momentum**](./man/MomentumMain.md), 
3) [**Mass** and **Compositon**](./man/AdvectMain.md). 

`GeoModBox.jl` includes a series of [exercises](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/) and [examples](./man/Examples.md) of geodynamically well-defined problems. The exercises are provided as Jupyter notebooks for students to complete. The theoretical background is documented here.

The solvers for each governing equation can be used separately or in combination for dimensional or non-dimensional problems, with only minimal modifications when calling the functions. Some typical initial conditions, such as a linearly increasing temperature, are predefined and can be called using [specific functions](./man/Ini.md).

## Staggered Finite Difference

To properly solve the governing equations, a staggered finite difference scheme is chosen for the *energy* and *momentum* equations. A staggered grid enables a correct and straightforward implementation of boundary conditions and ensures conservation of stress between nodes in cases of variable viscosity. This requires certain parameters to be defined on different grids.

Here, temperature, density, pressure, normal deviatoric stresses, and heat production rate are defined on the *centroids*. The deviatoric shear stresses are defined on the *vertices*, and velocities are defined between the *vertices*. Viscosity is required on both.

For further details on the implementation in `GeoModBox.jl`, see [here](./man/GESolution.md).

## Energy Conservation Equation

In geodynamics, the energy is described by the temperature and needs to be conserved within a closed system. Here, we solve the *temperature conservation equation*, or *temperature equation*, using an *operator splitting* method, that is, we first solve the *advective* part of the temperature equation, followed by the *diffusive* part. 

### [Heat Diffusion Equation](./man/DiffMain.md)

`GeoModBox.jl` provides several finite difference schemes for solving the *diffusive part* of the time-dependent or steady-state temperature equation, including radioactive heating, in both [1-D](./man/DiffOneD.md) and [2-D](./man/DiffTwoD.md). The solvers are located in [src/HeatEquation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/). Currently, only *Dirichlet* and *Neumann* thermal boundary conditions are supported. Most functions assume constant thermal parameters (with the exception of the 1-D solvers and the 2-D defect correction solver).

### [Heat Advection Equation](./man/AdvectMain.md)

```GeoModBox.jl``` provides various methods to advect properties within the model domain. The routines are structured so that any property defined on *centroids* (including *ghost nodes* at all boundaries) can be advected using the described solvers. Using passive tracers, one may choose to advect either the absolute temperature or the phase ID.

## [Momentum Conservation Equation](./man/MomentumMain.md)

On geological timescales, Earth's mantle and lithosphere deform slowly due to their high viscosity, allowing us to neglect inertial forces. This simplifies the Navier-Stokes equation into the **Stokes equation**. `GeoModBox.jl` provides two main methods to solve the Stokes equation in [1-D](./man/MomentumOneD.md) and [2-D](./man/MomentumTwoD.md): the direct method and the defect correction method, applicable for both constant and variable viscosity fields. Velocity and pressure are defined on a staggered grid, and ghost nodes are included to ensure proper implementation of free-slip and no-slip boundary conditions. 

## [Benchmarks and Examples](./man/Examples.md)

The following are visualizations of selected examples provided by `GeoModBox.jl`. For further details, refer to the documentation linked in each title.

### [Gaussian Temperature Diffusion](./man/examples/GaussianDiffusion2D.md)

![GaussianDiffusion](./assets/Gaussian_Diffusion_CNA_nx_100_ny_100.gif)

**Figure 1. Gaussian Diffusion.** Time-dependent, diffusive solution of a 2-D Gaussian temperature anomaly at a resolution of 100 × 100, using the [Crank-Nicholson approach](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/2Dsolvers.jl), compared to the analytical solution.  
Top Left: 2-D temperature field with numerical isotherms (solid black) and analytical isotherms (dashed yellow).  
Top Right: Total deviation from the analytical solution.  
Bottom Left: 1-D y-profile along $x = 0$.  
Bottom Right: Root Mean Square (RMS) total deviation over time.

![GDResTest](./assets/Gaussian_ResTest.png)

**Figure 2. Resolution test.** Maximum RMS error $\varepsilon$, maximum temperature, and mean temperature for various finite difference schemes and resolutions for the diffusion example shown above.

---

### [Rigid-Body-Rotation](./man/examples/Advection2D.md)

![RigidBodyI](./assets/2D_advection_circle_RigidBody_upwind_100_100_nth_1.gif)

![RigidBodyII](./assets/2D_advection_circle_RigidBody_semilag_100_100_nth_1.gif)

![RigidBodyIII](./assets/2D_advection_circle_RigidBody_tracers_100_100_nth_1.gif)

**Figure 3. Rigid-Body-Rotation.** Time-dependent advection of a rotating circular temperature anomaly using the **upwind (top)**, **semi-Lagrangian (middle)**, and **tracer (bottom)** methods on a 100 × 100 grid. Within a circular region, the velocity field follows rigid rotation; outside, it is zero. Temperature for tracers is interpolated to the grid for visualization but not updated on the tracers.

---

### [Falling Block](./man/examples/FallingBlockBenchmark.md)

![FallingBlockTD](./assets/Falling_block_ηr_0.0_tracers.gif)

**Figure 4. Isoviscous Falling Block.** Time-dependent simulation of an isoviscous falling block at 50 × 50 resolution with 9 tracers per cell. The solver handles variable viscosities. Tracers advect the phase ID, which is used to interpolate density and viscosity on centroids and vertices, respectively.

![FBSinkinVeloc](./assets/FallingBlock_SinkingVeloc_tracers.png)

**Figure 5. Falling Block Sinking Velocity.** Block sinking velocity vs. initial viscosity ratio $\eta_r$, using the same setup as above. 

![FBFinalStage](./assets/FallingBlock_FinalStage_tracers.png)

**Figure 6. Falling Block Benchmark.** Tracer distribution at the final stage for selected viscosity ratios $\eta_r \ge 0$.

--- 

### [Rayleigh-Taylor Instability](./man/examples/RTI.md)

![RTIani](./assets/RTI_ηr_-6.0_tracers_DC.gif)

**Figure 7. Rayleigh-Taylor Instability.** Evolution of two-layered Rayleigh-Taylor instability. 

---

### [Thermal Convection](./man/examples/MixedHeatedConvection.md)

![BHTC](./assets/Bottom_Heated_1.0e6_400_100_lineara.gif)

**Figure 8. Bottom-Heated, Isoviscous Convection for Ra = $10^6$, resolution 400 × 100.**  
TOP: Transient temperature field with velocity vectors.  
BOTTOM: Horizontally averaged temperature–depth profiles at each time step.  
Solvers: defect correction (momentum), semi-Lagrangian (advection), Crank-Nicolson (heat diffusion).  
Boundary conditions: Dirichlet (top/bottom), Neumann (sides), free-slip (velocity, all sides).

![IHTC](./assets/Internally_Heated_1.0e6_400_100_lineara.gif)

**Figure 9. Internally Heated Convection for $Ra_Q = 1.5 \cdot 10^6$, resolution 400 × 100.**  
Same setup as above, but with Neumann boundary at the bottom (zero heat flux) and constant internal volumetric heat production $Q \approx 15$.

![MHTC](./assets/Mixed_Heated_1.0e6_400_100_lineara.gif)

**Figure 10. Mixed-Heated Convection for Ra = $...$, resolution 400 × 100.**  
Combination of the above two setups (bottom heating + internal heating).

------------------

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://geosci-ffm.github.io/GeoModBox.jl/)
[![Unit Tests](https://github.com/GeoSci-FFM/GeoModBox.jl/actions/workflows/UnitTests.yml/badge.svg)](https://github.com/GeoSci-FFM/GeoModBox.jl/actions/workflows/UnitTests.yml)

# GeoModBox.jl

The **Geod**ynamic **Mod**elling Tool**Box** is a Julia package primarily intended for teaching purposes. It provides various finite difference, staggered discretization schemes to numerically solve the governing equations of two-dimensional geodynamic problems. These include the conservation equations of:

1) [**Energy**](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/DiffMain/), 
2) [**Momentum**](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/MomentumMain/), 
3) [**Mass** and **Composition**](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/AdvectMain/). 

`GeoModBox.jl` includes a series of [exercises](./exercises/) and [examples](./examples/) of geodynamically well-defined problems. The exercises are provided as Jupyter notebooks for students to complete. The theoretical background is documented here.

The solvers for each governing equation can be used separately or in combination for dimensional or non-dimensional problems, with only minimal modifications when calling the functions. For more information on how to use the individual functions please see the [list of functions](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/listoffunctions/) or individual [examples](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Examples/). Some typical initial conditions, such as a linearly increasing temperature, are predefined and can be called using [specific functions](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/Ini/). In the following a brief explanation is given regarding the governing equations and the numerical method to solve them within the `GeoModBox.jl`. For more detailed information see the individual documentations. 

## Installation

`GeoModBox.jl` can be installed directly through the Julia package manager or cloned from GitHub to access the complete repository, including the examples, exercises, and documentation.

> **Important:** `GeoModBox.jl` uses separate Julia environments for the core package, examples, exercises, and documentation. When running examples or exercises, the corresponding `examples` or `exercises` environment must be activated first.

For detailed instructions on installing Julia and `GeoModBox.jl`, setting up the required environments, and using the package as a user, student, or developer, see the [Installation](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/Installation/) guide.

## Try GeoModBox.jl online

### Google Colab

Selected `GeoModBox.jl` exercises and examples can be executed directly in Google Colab without a local Julia installation.

- Exercise 02 – 1-D Heat Diffusion (Explicit)</br>
    - Student exercise [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeoSci-FFM/GeoModBox.jl/blob/main/colab/exercises/master/02_1D_Heat_explicit_en.ipynb)
    - Master solution [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeoSci-FFM/GeoModBox.jl/blob/main/colab/exercises/students/02_1D_Heat_explicit_en.ipynb)

- Exercise 10 – 2-D Falling Block (time-dep.)</br>
    - Student exercise [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeoSci-FFM/GeoModBox.jl/blob/main/colab/exercises/master/10_2D_Falling_Block_td_en.ipynb)
    - Master solution [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeoSci-FFM/GeoModBox.jl/blob/main/colab/exercises/students/10_2D_Falling_Block_td_en.ipynb)

- Bottom heated, isoviscous thermal Convection [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GeoSci-FFM/GeoModBox.jl/blob/main/colab/examples/BottomHeated.ipynb)

## Staggered Finite Difference

To properly solve the governing equations, a staggered finite difference scheme is chosen for the *energy* and *momentum* equations. A staggered grid enables a correct and straightforward implementation of boundary conditions and ensures conservation of stress between nodes in cases of variable viscosity. This requires certain parameters to be defined on different grids. For more information regarding the physical and numerical background, please refer to [this](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/GESolution/).

Within the `GeoModBox.jl`, temperature, density, pressure, normal deviatoric stresses, and heat production rate are defined on the *centroids*. The deviatoric shear stresses are defined on the *vertices*, and velocities are defined between the *vertices*. Viscosity is required on both.

For further details on the implementation in `GeoModBox.jl`, see the individual documentations for each governing equation. 

## Energy Conservation Equation

In geodynamics, the energy is described by the temperature and needs to be conserved within a closed system. Within the `GeoModBox.jl`, the *temperature conservation equation*, or *temperature equation*, is solved using an *operator splitting* method, that is, first the *advective* part of the temperature equation is solved, followed by the *diffusive* part. 

### [Heat Diffusion Equation](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/DiffMain/)

`GeoModBox.jl` provides several finite difference schemes for solving the *diffusive part* of the time-dependent or steady-state temperature equation, including radioactive heating, in both [1-D](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/DiffOneD/) and [2-D](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/DiffTwoD/). The solvers are located in [src/HeatEquation](./src/HeatEquation/). Currently, only *Dirichlet* and *Neumann* thermal boundary conditions are supported. Most functions assume constant thermal parameters (with the exception of the 1-D solvers and the 2-D, iterative implicit solver, called **iterative defect correction method**).

### [Heat Advection Equation](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/AdvectMain/)

`GeoModBox.jl` provides various methods to advect properties within the model domain. The routines are structured so that any property defined on *centroids* (including *ghost nodes* at all boundaries) can be advected using the described solvers. Using passive tracers, one may choose to advect either the absolute temperature or the phase ID.

## [Momentum Conservation Equation](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/MomentumMain/)

On geological timescales, Earth's mantle and lithosphere deform slowly due to their high viscosity, allowing us to neglect inertial forces. This simplifies the Navier-Stokes equation into the **Stokes equation**. `GeoModBox.jl` provides two main methods to solve the Stokes equation in [1-D](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/MomentumOneD/) and [2-D](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/MomentumTwoD/): the direct method and the defect correction method, applicable for both constant and variable viscosity fields. Velocity and pressure are defined on a staggered grid, and ghost nodes are included to ensure proper implementation of free-slip and no-slip boundary conditions.

## [Benchmarks and Examples](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Examples/)

The following are visualizations of selected examples provided by `GeoModBox.jl`. For further details, refer to the documentation linked in each title.  Each runtime can be found in the documentation of the [examples](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Examples/) and [exercises](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/exercises/Exercises/).

### [Gaussian Temperature Diffusion](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Diffusion/GaussianDiffusion2D/)

<img src="./docs/src/assets/examples/Diffusion/Gaussian_Diffusion_CN_nx_120_ny_120.gif" alt="drawing" width="600"/>

**Figure 1. Gaussian Diffusion.** Time-dependent, diffusive solution of a 2-D Gaussian temperature anomaly at a resolution of 120 × 120, using the special solver with the [Crank-Nicholson discretization](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/theory/DiffMain/), compared to the analytical solution. a) 2-D temperature field with numerical isotherms (solid black) and analytical isotherms (dashed yellow). b) Total deviation from the analytical solution. c) 1-D y-profile along $x = 0$. d) Root Mean Square (RMS) total deviation over time.

<img src="./docs/src/assets/examples/Diffusion/Gaussian_ResTest.png" alt="drawing" width="600"/>

**Figure 2. Resolution test.** a) Maximum RMS error $\varepsilon$, b) maximum temperature, and c) mean temperature for various finite difference schemes and resolutions using the special case, linear solver (single left-matrix division) for the diffusion example shown above.

---

### [Rigid-Body-Rotation](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Advection/Advection2D/)

<img src="./docs/src/assets/examples/Advection/2D_advection_circle_RigidBody_upwind_100_100_nth_1.gif" alt="drawing" width="600"/> 

<img src="./docs/src/assets/examples/Advection/2D_advection_circle_RigidBody_semilag_100_100_nth_1.gif" alt="drawing" width="600"/>

<img src="./docs/src/assets/examples/Advection/2D_advection_circle_RigidBody_markers_100_100_nth_1.gif" alt="drawing" width="600"/> 
<br>

**Figure 3. Rigid-Body-Rotation.** Time-dependent advection of a rotating circular temperature anomaly using the **upwind (top)**, **semi-Lagrangian (middle)**, and **tracer (bottom)** methods on a 100 × 100 grid. Within a circular region, the velocity field follows rigid rotation; outside, it is zero.

---

### [Falling Block](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Stokes/FallingBlockBenchmark/)

<img src="./docs/src/assets/examples/Stokes/Falling_block_ηr_0.0_tracers_DC.gif" alt="drawing" width="600"/>

**Figure 4. Isoviscous Falling Block.** Time-dependent simulation of an isoviscous falling block at 50 × 50 resolution with 9 tracers per cell. The solver handles variable viscosities. Tracers advect the phase ID, which is used to interpolate density and viscosity on centroids and vertices, respectively.

<img src="./docs/src/assets/examples/Stokes/FallingBlock_SinkingVeloc_tracers_direct_arithmetic.png" alt="drawing" width="600"/>

**Figure 5. Falling Block Sinking Velocity.** Block sinking velocity vs. initial viscosity ratio $\eta_r$, using the same setup as above. 

<img src="./docs/src/assets/examples/Stokes/FallingBlock_FinalStage_tracers_direct_arithmetic.png" alt="drawing" width="600"/>

**Figure 6. Falling Block Benchmark.** Tracer distribution at the final stage for selected viscosity ratios $\eta_r \ge 0$.

--- 

### [Rayleigh-Taylor Instability](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Stokes/RTI/)

#### [RTI Growth-rate Benchmark (Ramberg, 1968)](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Stokes/RTI_growth_rate/)

<img src="./docs/src/assets/examples/Stokes/RTI_ηr_-6.0_tracers_DC_MarkerInterpolation_arith.gif" alt="drawing" width="600"/>

**Figure 7. Rayleigh-Taylor Instability.** Evolution of two-layered Rayleigh-Taylor instability. 

<img src="./docs/src/assets/examples/Stokes/RTI_Growth_Rate_nmx_5_nmy_5_MarkerInterpolation_arith.png" width="600">

**Figure 8.** Growth rate of an initial cosinusoidal perturbation in a two-layer system across various wavelengths $\lambda$. The growth rate is arbitrarily scaled using $b_1$ and $b_2$ for visualization, following the approach of Gerya (2019). The lines are the analytical solutions for different viscosity ratios $\eta_r$ and the markers show the corresponding numerical results for models with decreasing amplitudes (black - 100 m, red - 10 m, yellow - 1 m).

#### [RTI Benchmark (Van Keken et al., 1997)](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Stokes/VanKekenBenchmark/)

<img src="./docs/src/assets/examples/Stokes/VanKeKen_Benchmark_ηr_0.0_tracers_DC_arith.gif" width="600">

**Figure 9.** Evolution of the dimensional Van Keken benchmark. The panels show the density field, tracer distribution, viscosity, and velocity magnitude together with the velocity vectors. The initially perturbed interface evolves into the characteristic Rayleigh-Taylor instability as the denser lower material sinks beneath the lighter upper layer.

---

### [Viscous Inclusion Benchmark](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Stokes/ViscousInclusion/)

<img src="./docs/src/assets/examples/Stokes/ViscousInclusion_Summary.png" width="600">

**Figure 10.** Comparison between the numerical and analytical solutions for the horizontal velocity, vertical velocity, and pressure fields. The third column shows the corresponding relative error distributions. The largest errors are localized at the viscosity interface, where the material properties are discontinuous, while excellent agreement is obtained throughout the remainder of the domain.

---

### Thermal Convection Models

#### [Blankenbach Benchmark (Blankenbach et al., 1989)](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/exercises/13_Blankenbach_Benchmark/)

<img src="./docs/src/assets/exercises/13a.gif" width="700">

**Figure 11.** Isoviscous, bottom-heated thermal convection for $Ra_b = 10^6$ with a resolution of 100×100.  
The initial condition is a sinusoidally perturbed conductive temperature field.  
The background color shows the non-dimensional temperature, overlaid by temperature isolines (every 0.05) and centroid velocity vectors. Heat diffusion is solved using the defect correction with a **Crank–Nicolson** discretization, the Stokes equation using the **defect correction** method, and temperature advection with the **semi-Lagrangian** method.  

<img src="./docs/src/assets/exercises/13g.png" width="700">

**Figure 12.** Summary of the resolution study for a basal Rayleigh number of $Ra_b = 10^6$. (a) Final dimensionless temperature field with superimposed velocity vectors and temperature contours for the third model. (b) Vertical temperature profile through the center of the model. The black squares denote the published locations of the minimum and maximum temperatures reported by Blankenbach et al. (1989). (c) Evolution of the Nusselt number ($Nu$) and root-mean-square velocity ($V_{\mathrm{RMS}}$) toward steady state. The dashed red lines indicate the benchmark values. (d) Grid convergence of the Nusselt number, root-mean-square velocity, and mean temperature as a function of the inverse number of grid cells, $1/(n_x n_y)$. The dashed red lines correspond to the benchmark values.

<img src="./docs/src/assets/examples/Convection/Blankenbach_VarEta_1.00e+04_100_100_blankenbach.gif" width="700">

**Figure 13.** Evolution of the two-dimensional Blankenbach benchmark with temperature-dependent viscosity for a resolution of $100 \times 100$. The left panel shows the dimensionless temperature field with superimposed velocity vectors, while the right panel displays the logarithm of the normalized viscosity, $\log_{10}(\eta/\eta_0)$. As the initially conductive state becomes unstable, cold, highly viscous material accumulates beneath the upper boundary, whereas the hot interior develops a pronounced low-viscosity region that promotes vigorous convection.

#### [Mixed Thermal Convection Models](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Convection/Overview_Convection/)

##### [Bottom Heated Convection](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Convection/BottomHeatedConvection.md)

<img src="./docs/src/assets/examples/Convection/Bottom_Heated_1.0e6_400_100_lineara.gif" width="700">

**Figure 14. Bottom-Heated, Isoviscous Convection for $Ra_b = 10^6$, resolution 400 × 100.**  
Solvers: defect correction(momentum), semi-Lagrangian (advection), defect correction + Crank-Nicolson discretization (heat diffusion).  
Boundary conditions: Dirichlet (top/bottom), Neumann (sides), free-slip (velocity, all sides).

##### [Internally Heated Convection](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Convection/InternallyHeatedConvection/)

<img src="./docs/src/assets/examples/Convection/Internally_Heated_1.0e6_400_100_lineara.gif" width="700">

**Figure 15. Internally Heated Convection for $Ra_Q \approx 1.5 \cdot 10^7$, resolution 400 × 100.**  
Same setup as above, but with Neumann boundary at the bottom (zero heat flux) and constant internal volumetric heat production $Q \approx 15$.

##### [Mixed Heated Convection](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Convection/MixedHeatedConvection/)

<img src="./docs/src/assets/examples/Convection/Mixed_Heated_1.0e6_400_100_lineara.gif" width="700">

**Figure 16. Mixed-Heated Convection for a basal and internal-heating Rayleigh number of $Ra_b = 10^6, Ra_Q \approx 1.5 \cdot 10^7$, resolution 400 × 100.**  
Combination of the above two setups (bottom heating + internal heating).

#### [Bottom Heated with temperature-dependent viscosity](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/Convection/BottomHeatedConvection_VE/)

<img src="./docs/src/assets/examples/Convection/Bottom_Heated_VE_1.0e6_400_100_lineara.gif" width="700">

**Figure 17. Bottom-Heated, variable viscosity convection for $Ra_b = 10^6$, resolution 400 × 100.** Same setup as above, but with a temperature-dependent, Arrhenius-like viscosity resulting in a viscosity contrast of five orders of magnitude in a non-dimensional temperature range of 0-1.

---

### [Thermo-mechanical Shear Localization](https://geosci-ffm.github.io/GeoModBox.jl/dev/man/examples/StrainLocalization/)

<img src="./docs/src/assets/examples/ShearHeating/ShearHeatingBands_2D.gif" width="700">

**Figure 18.** Evolution of thermo-mechanical shear localization in the pure-shear benchmark using shear heating. The model employs the general defect-correction energy solver with a Crank–Nicolson discretization ($\theta = 0.5$), tracer-based temperature advection, and a grid resolution of $200 \times 100$ cells. Temperature-dependent, non-linear dislocation-creep viscosities are combined using arithmetic averaging for the phase mixing, centroid-to-vertex interpolation, and marker-to-grid interpolation. The animation shows the evolution of the second invariant of the strain-rate field (background colors), the weak inclusion (black contour), the evolving shear band (white contour), and the fixed diagnostic profile (red line) used to quantify strain-rate amplification, shear-band thickness, and shear-heating-induced temperature increase.

<img src="./docs/src/assets/examples/ShearHeating/FinalBenchmarkFigure_fixed_arithmetic_arithmetic.png" width="700">

**Figure 19.** Evolution and quantitative characterization of thermo-mechanical shear localization for the fixed-profile diagnostic with arithmetic phase and vertex viscosity averaging. Panels (a)–(c) show the second invariant of the strain-rate field at approximately 5%, 15%, and 25% bulk shortening, including the weak inclusion and the profile used to evaluate localization. Panels (d)–(f) show the the strain-rate amplification, temperature increase, shear-band orientation, and shear-band thickness as functions of bulk shortening. Line color distinguishes the upwind, semi-Lagrangian, and tracer advection methods, while line style distinguishes backward Euler (θ = 0), Crank–Nicolson (θ = 0.5), and forward Euler (θ = 1) within the general defect-correction formulation.

------------------

<!-- # References

Gerya, T. (2019). Introduction to numerical geodynamic modelling. Cambridge University Press.

Spiegelman, M. (2004). Myths and methods in modeling. Columbia University Course Lecture Notes, available online at http://www. ldeo. columbia. edu/~ mspieg/mmm/course. pdf, accessed, 6, 2006.

Becker, T.W., and Kaus, B.J.P., 2016, Numerical modeling of Earth systems, an introduction to computational methods with focus on solid Earth applications of continuum mechanics: University of Southern California Lecture Notes, http://www-udc.ig.utexas.edu/external/becker/preprints/Geodynamics557.pdf (last accessed June 2025).

W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes 1986, (Cambridge Univ. Press, Cambridge, 1986).

Becker, T., & Faccenna, C. (2025). Tectonic Geodynamics. Princeton University Press.

# References

Spiegelman, M. (2004). Myths and methods in modeling. Columbia University Course Lecture Notes, available online at http://www. ldeo. columbia. edu/~ mspieg/mmm/course. pdf, accessed, 6, 2006.

W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes 1986, (Cambridge Univ. Press, Cambridge, 1986).
-->

<!--
- Blanckenbach
- Channel Flow (2D)
- Falling Block, check! 
- Gauss Diffusion, check! 
- RTI 
- Rigid Body Rotation, check! 
- Viscous Inclusion
-->

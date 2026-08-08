# Exercises

`GeoModBox.jl` contains 13 interactive exercises designed to introduce and practice the numerical solution of partial differential equations. The focus is on the **finite difference method** for solving governing equations in geodynamics.

The sequence of exercises begins with a simple 1-D Euler advection problem and proceeds to the discretization and solution of the 1-D heat diffusion equation using different numerical schemes. The heat equation is then extended to two dimensions, allowing students to analyze and compare the advantages and disadvantages of various approaches.

Next, the pure 1-D advection equation is revisited and solved with different numerical schemes, followed by coupling of the advection and diffusion equations using the operator-splitting method.

The exercises then move on to the 1-D Stokes equation, solved using both direct and iterative defect correction methods, before progressing to the steady-state and time-dependent solutions of a full Stokes problem.

The final section addresses the three fundamental conservation laws in geodynamics:  
1. Conservation of energy,  
2. Conservation of momentum, and  
3. Conservation of mass.  

These exercises investigate 2-D thermal convection, first in dimensional form and then in a scaled (non-dimensional) version.  

The series concludes with a reproduction of the **Blankenbach benchmark**.

Below, the runtimes for each exercise and their individual options are listed. All exercises were performed on a single CPU: *AMD Ryzen 7 7735U with Radeon Graphics*.

| Exercise | Type | Total Runtime [s] |
|:----------|:------|:----------------|
| 1 | [Euler Advection](./01_Euler_Advection.md) | 5.741 |
| 2 | [1-D Heat Diffusion (explicit)](./02_1D_Heat_explicit.md) | 1) 1.82 s |
|   |   | 2) 1.75 s  |
|   |   | 3) 1.86 s |
| 3 | [1-D Heat Diffusion (implicit)](./03_1D_Heat_implicit.md) | 1) 1.48 s  |
|   |   | 2) 1.57 s  |
|   |   | 3) 2.55 s |
| 4 | [2-D Steady State Heat Equation](./04_2D_Diffusion_Stationary.md) | 3.69 s |
| 5a | [2-D Time-dep. Heat Equation (Plume)](./05_2D_Diffusion_TD_Plume.md) | 1) Special Solver (Explicit): 19.74 s |
|    |   | 2) Special Solver (Implicit): 55.29 s |
|    |   | 3) General Solver (Explicit): 24.75 s |
|    |   | 4) General Solver (Implicit): 58.44 s |
|    |   | 5) General Solver (CNA): 59.35 s |
| 5b | [2-D Time-dep. Heat Equation (Sill)](./05_2D_Diffusion_TD_Sill.md) | 1) Special Solver (Explicit): 17.08 s |
|    |   | 2) Special Solver (Implicit): 56.61 s |
|    |   | 3) General Solver (Explicit): 31.24 |
|    |   | 4) General Solver (Implicit): 60.42 s |
|    |   | 5) General Solver (CNA): 58.9 s |
| 6 | [1-D Advection Scheme](./06_1D_Advection.md) | 1) FTCS: 4.866 |
|   |   | 2) Upwind: 5.047 |
|   |   | 3) Downwind: 4.68 |
|   |   | 4) LAX: 5.202 |
|   |   | 5) SLF: 5.241 |
|   |   | 6) Semilag: 9.886 |
|   |   | 7) Tracers: 9.566 |
| 7 | [2-D Energy Equation](./07_2D_Energy_Equation.md) | 1) Upwind+explicit: 9.11 s |
|   |   | 2) Upwind+implicit: 19.04 s |
|   |   | 3) Upwind+CNA: 16.35 s |
|   |   | 4) Upwind+ADI: 144.68 s |
|   |   | 5) Upwind+dc: 23.23 s |
|   |   | 6) SL+explicit: 19.1 s |
|   |   | 7) SL+implicit: 17.11 s |
|   |   | 8) SL+CNA: 16.55 s |
|   |   | 9) SL+ADI: 148.18 s |
|   |   | 10) SL+dc: 18.61 s |
|   |   | 11) Tracers+explicit: 120.41 s |
|   |   | 12) Tracers+implicit: 143.14 s |
|   |   | 13) Tracers+CNA: 143.1 s |
|   |   | 14) Tracers+ADI: 246.73 s |
|   |   | 15) Tracers+dc: 150.38 s |
| 8 | [1-D Stokes Equation](./08_1D_Stokes.md) | 671 ms |
| 9 | [2-D Falling Block (steady state)](./09_2D_Falling_Block.md) | 1.22 s |
| 10 | [2-D Falling Block (time-dep.)](./10_2D_Falling_Block_td.md) | 1) Upwind: 4.89 s |
|    |   | 2) SLF: 5.31 s |
|    |   | 3) SL: 4.94 s |
|    |   | 4) Tracers: 7.96 s |
| 11¹ | [Thermal Convection (dim)](./11_2D_Thermal_Convection.md), **Resolution: 150×50** | **Ra = 1e4** (Diff+Adv+Momentum) |
|     |   | 1) Explicit+upwind+direct: 353.2 s |
|     |   | 2) Implicit+upwind+direct: 438.67 s |
|     |   | 3) CNA+upwind+direct: 440.02 s |
|     |   | 4) ADI+upwind+direct: 966.84 s |
|     |   | 5) DC+upwind+direct: 383.1 s |
|     |   | 6) Explicit+semilag+direct: 409.16 s | 
|     |   | 7) Implicit+semilag+direct: 432.15 s | 
|     |   | 8) CNA+semilag+direct: 420.52 s |
|     |   | 9) ADI+semilag+direct: 976.48 s | 
|     |   | 10) DC+semilag+direct: 456.91 s | 
|     |   | 11) Explicit+upwind+DC: 384.52 s |
|     |   | 12) Implicit+upwind+DC: 454.44 s |
|     |   | 13) CNA+upwind+DC: 456.66 s |
|     |   | 14) ADI+upwind+DC: 1007.67 s |
|     |   | 15) DC+upwind+DC:  412.36 s|
|     |   | 16) Explicit+semilag+DC: 428.12 s |
|     |   | 17) Implicit+semilag+DC: 461.63 s |
|     |   | 18) CNA+semilag+DC: 453.25 s |
|     |   | 19) ADI+semilag+DC: 1045.37 s |
|     |   | 20) DC+semilag+DC: 471.471 s |
|     |   | **Ra = 1e5** (Diff+Adv+Momentum) |
|     |   | 21) CNA+upwind+dc: 367.553 s |
|     |   | 22) Explicit+semilag+dc: 352.524 s |
|     |   | 23) CNA+semilag+dc: 357.559 s |
|     |   | 24) dc+semilag+dc: 358.312 |
|     |   | **Ra = 1e6** (Diff+Adv+Momentum)** |
|     |   | 25) CNA+upwind+dc: 1299.33 s |
|     |   | 26) Explicit+semilag+dc: 1205.31 s |
|     |   | 27) CNA+semilag+dc: 1288.58 s |
|     |   | 28) dc+semilag+dc: 1311.65 s|
| 12² | [Thermal Convection (scaled)](./12_2D_Thermal_Convection_scaled.md), **Resolution: 150×50** | **Ra = 1e4** (Diff+Adv+Momentum) |
|     |   | 1) Explicit+upwind+direct: 149.4 s |
|     |   | 2) Implicit+upwind+direct: 166.21 s |
|     |   | 3) CNA+upwind+direct: 167.71 s |
|     |   | 4) ADI+upwind+direct: 602.02 s |
|     |   | 5) DC+upwind+direct: 172.58 s |
|     |   | 6) Explicit+semilag+direct: 149.8 s | 
|     |   | 7) Implicit+semilag+direct: 174.48 s | 
|     |   | 8) CNA+semilag+direct: 175.7 s |
|     |   | 9) ADI+semilag+direct: 662.41 s | 
|     |   | 10) DC+semilag+direct: 177.45 s | 
|     |   | 11) Explicit+upwind+DC: 150.88 s |
|     |   | 12) Implicit+upwind+DC: 169.41 s |
|     |   | 13) CNA+upwind+DC: 173.95 s |
|     |   | 14) ADI+upwind+DC: 606.75 s |
|     |   | 15) DC+upwind+DC: 172.08 s |
|     |   | 16) Explicit+semilag+DC: 155.17 s |
|     |   | 17) Implicit+semilag+DC: 178.08 s |
|     |   | 18) CNA+semilag+DC: 178.45 s |
|     |   | 19) ADI+semilag+DC: 659.53 s |
|     |   | 20) DC+semilag+DC: 184.13 s |
|     |   | **Ra = 1e5** (Diff+Adv+Momentum) |
|     |   | 21) CNA+upwind+dc: 137.31 s |
|     |   | 22) Explicit+semilag+dc: 117.15  |
|     |   | 23) CNA+semilag+dc: 137.44 s |
|     |   | 24) dc+semilag+dc: 159.57 s |
|     |   | **Ra = 1e6** (Diff+Adv+Momentum)** |
|     |   | 25) CNA+upwind+dc: 365.57 s |
|     |   | 26) Explicit+semilag+dc: 302.45 s |
|     |   | 27) CNA+semilag+dc: 383.05 s |
|     |   | 28) dc+semilag+dc: 384.23 s|
| 13 | [Blankenbach Benchmark](./13_Blankenbach_Benchmark.md) | **Resolution: 50×50** |
|    |   | Ra = 1e4: 91.6504 s |
|    |   | Ra = 1e5: 91.1549 s |
|    |   | Ra = 1e6: 321.539 s |
|    |   | **Resolution: 100×100** |
|    |   | Ra = 1e4:  563.994 s |
|    |   | Ra = 1e5:  297.469 s |
|    |   | Ra = 1e6:  1069.46 s |
|    |   | **Resolution Tests:** |
|    |   | Ra = 1e4: 160.83 s | 
|    |   | Ra = 1e5: 1142.15 s | 
|    |   | Ra = 1e6: 36283.8 s | 

¹: For the higher $Ra$ number models only the fastest and most accurate (see [here](../examples/Diffusion/GaussianDiffusion2D.md)) combinations of solvers from the low $Ra$ case haven been run. The runtime was determined by executing the code in a seperate Julia script. 

²: Besides the non-dimensionalization the exercise was optimized numerically to reduce the runtime. The runtime was determined by executing the code in a seperate Julia script. 
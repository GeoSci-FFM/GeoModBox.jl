# Exercises

`GeoModBox.jl` contains 13 interactive exercises designed to introduce and practice the numerical solution of partial differential equations. The focus is on the **finite difference method** for solving governing equations in geodynamics.

The sequence of exercises begins with a simple 1-D Euler advection problem and proceeds to the discretization and solution of the 1-D heat diffusion equation using different numerical schemes. The heat equation is then extended to two dimensions, allowing students to analyze and compare the advantages and disadvantages of various approaches.

Next, the pure 1-D advection equation is revisited and solved with different numerical schemes, followed by coupling of the advection and diffusion equations using the operator-splitting method.

The exercises then move on to the 1-D Stokes equation, solved using both direct and iterative defect-correction methods, before progressing to the steady-state and time-dependent solutions of a full Stokes problem.

The final section addresses the three fundamental conservation laws in geodynamics:  
1. Conservation of energy,  
2. Conservation of momentum, and  
3. Conservation of mass.  

These exercises investigate 2-D thermal convection, first in dimensional form and then in a scaled (non-dimensional) version.  

The series concludes with a reproduction of the **Blankenbach benchmark**.

---

Below, the runtimes for each exercise and their individual options are listed.  
All exercises were performed on a single CPU: *AMD Ryzen 7 7735U with Radeon Graphics*.  

- The **first runtime** corresponds to the initial execution, including compilation and memory allocation.  
- The **runtime in parentheses** corresponds to a second execution, where compilation and allocation overhead are avoided.  

Test Table I: 

<table>
  <tr>
    <th>#</th>
    <th>Type</th>
    <th>Runtime [s]</th>
  </tr>
  <tr>
    <td>1</td>
    <td>Euler Advection</td>
    <td>5.741</td>
  </tr>
  <tr>
    <td>2</td>
    <td>1-D Heat Diffusion (explicit)</td>
    <td>
      1) 9.342 (1.14)<br>
      2) 9.975 (1.087)<br>
      3) 9.616 (1.017)
    </td>
  </tr>
  <tr>
    <td>3</td>
    <td>1-D Heat Diffusion (implicit)</td>
    <td>
      1) 12.641 (1.749)<br>
      2) 13.537 (1.720)<br>
      3) 12.157 (1.731)
    </td>
  </tr>
</table>

Test Table II: 

| Exercise | Type | Total Runtime [s] |
|:----------|:------|:----------------|
| 1 | Euler Advection | 5.741 |
| 2 | 1-D Heat Diffusion (explicit) | 1) 9.342 (1.14) |
|   |   | 2) 9.975 (1.087) |
|   |   | 3) 9.616 (1.017) |
| 3 | 1-D Heat Diffusion (implicit) | 1) 12.641 (1.749) |
|   |   | 2) 13.537 (1.720) |
|   |   | 3) 12.157 (1.731) |
| 4 | 2-D Steady State Heat Equation | 13.471 (2.069) |
| 5a | 2-D Time-dep. Heat Equation (Plume) | 1) Explicit: 18.632 (4.799) |
|    |   | 2) Implicit: 78.201 (60.828) |
| 5b | 2-D Time-dep. Heat Equation (Sill) | 1) Explicit: 17.856 (4.566) |
|    |   | 2) Implicit: 39.626 (22.679) |
| 6 | 1-D Advection Scheme | 1) FTCS: 11.644 (4.019) |
|   |   | 2) Upwind: 11.03 (4.042) |
|   |   | 3) Downwind: 11.465 (3.343) |
|   |   | 4) LAX: 11.272 (3.405) |
|   |   | 5) SLF: 11.312 (3.48) |
|   |   | 6) Semilag: 14.318 (3.604) |
|   |   | 7) Tracers: 15.662 (4.944) |
| 7 | 2-D Energy Equation | 1) Upwind+explicit: 17.468 (2.979) |
|   |   | 2) Upwind+implicit: 31.441 (10.842) |
|   |   | 3) Upwind+CNA: 32.065 (12.222) |
|   |   | 4) Upwind+ADI: 121.665 (99.807) |
|   |   | 5) Upwind+dc: 30.45 (11.401) |
|   |   | 6) SLF+explicit: 22.002 (2.872) |
|   |   | 7) SLF+implicit: 28.751 (8.911) |
|   |   | 8) SLF+CNA: 29.876 (9.28) |
|   |   | 9) SLF+ADI: 120.018 (98.74) |
|   |   | 10) SLF+dc: 28.179 (9.379) |
|   |   | 11) SL+explicit: 22.969 (6.334) |
|   |   | 12) SL+implicit: 30.634 (9.682) |
|   |   | 13) SL+CNA: 32.415 (11.566) |
|   |   | 14) SL+ADI: 133.871 (110.704) |
|   |   | 15) SL+dc: 35.364 (12.301) |
|   |   | 16) Tracers+explicit: 89.499 (65.913) |
|   |   | 17) Tracers+implicit: 108.464 (84.209) |
|   |   | 18) Tracers+CNA: 112.171 (89.219) |
|   |   | 19) Tracers+ADI: 186.71 (161.611) |
|   |   | 20) Tracers+dc: 115.737 (89.089) |
| 8 | 1-D Stokes Equation | 8.489 (0.212) |
| 9 | 2-D Falling Block (steady state) | 20.764 (0.891) |
| 10 | 2-D Falling Block (time-dep.) | 1) Upwind: 33.295 (3.68) |
|    |   | 2) SLF: 32.333 (3.629) |
|    |   | 3) SL: 37.537 (3.755) |
|    |   | 4) Tracers: 36.291 (6.12) |
| 11¹ | Thermal Convection (dim), **Resolution: 150×50** | **Ra = 1e4** (Diff+Adv+Momentum) |
|     |   | 1) Explicit+upwind+direct: 262.912 (1703) |
|     |   | 2) Implicit+upwind+direct: 283.031 (1709) |
|     |   | 3) CNA+upwind+direct: 284.944 (1708) |
|     |   | 4) ADI+upwind+direct: 645.181 (1708) |
|     |   | 5) DC+upwind+direct: 294.1 (1710) |
|     |   | 6) Explicit+SLF+direct: 354.35 (2341) |
|     |   | 7) Implicit+SLF+direct: 382.307 (2345) | 
|     |   | 8) CNA+SLF+direct: 386.764 (2344) |
|     |   | 9) ADI+SLF+direct: 892.411 (2344) | 
|     |   | 10) DC+SLF+direct: 492.168 (2346) | 
|     |   | 11) Explicit+semilag+direct: 334.941 (1743) | 
|     |   | 12) Implicit+semilag+direct: 367.501 (1746) | 
|     |   | 13) CNA+semilag+direct: 371.289 (1744) |
|     |   | 14) ADI+semilag+direct: 704.852 (1744) | 
|     |   | 15) DC+semilag+direct: 390.767 (1749) | 
|     |   | 16) Explicit+upwind+DC: 371.091 (1703) |
|     |   | 17) Implicit+upwind+DC: 358.323 (1709) |
|     |   | 18) CNA+upwind+DC: 338.862 (1708) |
|     |   | 19) ADI+upwind+DC: 654.069 (1708) |
|     |   | 20) DC+upwind+DC: 316.456 (1710) |
|     |   | 21) Explicit+SLF+DC: 394.471 (2341) |
|     |   | 22) Implicit+SLF+DC: 417.712 (2345) |
|     |   | 23) CNA+SLF+DC: 412.717 (2344) |
|     |   | 24) ADI+slf+DC: 890.714 (2344) |
|     |   | 25) DC+SLF+DC: 403.217 (2346) |
|     |   | 26) Explicit+semilag+DC: 286.046 (1743) |
|     |   | 27) Implicit+semilag+DC: 295.961 (1746) |
|     |   | 28) CNA+semilag+DC: 288.067 (1744) |
|     |   | 29) ADI+semilag+DC: 671.159 (1744) |
|     |   | 30) DC+semilag+DC: 357.894 (1749) |
|     |   | **Ra = 1e5** (Diff+Adv+Momentum) |
|     |   | 1) Explicit+semilag+dc: 489.869 (2433) |
|     |   | 2) CNA+semilag+dc: 515.539 (2448) |
|     |   | 3) CNA+upwind+dc: 1327.856 (8000) |
|     |   | **Ra = 1e6** (Diff+Adv+Momentum)** |
|     |   | 4) Explicit+semilag+dc: 1280.077 (8000) |
|     |   | 5) CNA+semilag+dc: 1324.29 (8000) |
|     |   | 6) CNA+upwind+dc: 1272.997 (8000) |
| 12 | Thermal Convection (scaled), **Resolution: 150×50** | **Ra = 1e4** (Diff+Adv+Momentum) |
|    |   | 1) Explicit+upwind+direct: 310.208 (1950) |
|    |   | 2) Implicit+upwind+direct: 333.82 (1957) |
|    |   | 3) CNA+upwind+direct: 331.791 (1955) |
|    |   | 4) ADI+upwind+direct: 739.821 (1955) |
|    |   | 5) DC+upwind+direct: 352.163 (1953) |
|    |   | 6) Explicit+SLF+direct: 437.435 (2704) |
|    |   | 7) Implicit+SLF+direct: 468.642 (2708) |
|    |   | 8) CNA+SLF+direct: 464.97 (2707) |
|    |   | 9) ADI+SLF+direct: 1013.196 (2707) |
|    |   | 10) DC+SLF+direct: 476.063 (2706) |
|    |   | 11) Explicit+semilag+direct: 358.84 (1999) |
|    |   | 12) Implicit+semilag+direct: 341.811 (2002) |
|    |   | 13) CNA+semilag+direct: 384.29 (2000) |
|    |   | 14) ADI+semilag+direct: 754.453 (2000) |
|    |   | 15) DC+semilag+direct: 359.703 (2001) |
|    |   | 16) Explicit+upwind+DC: 683.465 (1950) |
|    |   | 17) Implicit+upwind+DC: 614.899 (1957) |
|    |   | 18) CNA+upwind+DC:  649.386 (1955) |
|    |   | 19) ADI+upwind+DC: 1048.092 (1955) |
|    |   | 20) DC+upwind+DC:  640.199 (1953) |
|    |   | 21) Explicit+SLF+DC:  835.568 (2704) |
|    |   | 22) Implicit+SLF+DC:  862.991 (2708) |
|    |   | 23) CNA+SLF+DC:  868.07 (2707) |
|    |   | 24) ADI+slf+DC:  1402.997 (2707) |
|    |   | 25) DC+SLF+DC: 876.505 (2706) |
|    |   | 26) Explicit+semilag+DC: 613.413 (1999) |
|    |   | 27) Implicit+semilag+DC: 652.126 (2002) |
|    |   | 28) CNA+semilag+DC: 638.264 (2000) |
|    |   | 29) ADI+semilag+DC:  1027.153 (2000) |
|    |   | 30) DC+semilag+DC: 620. 945 (2001) <br/> <br/> |
|    |   | **Ra = 1e5** (Diff+Adv+Momentum)** |
|    |   | 31) Explicit+semilag+dc: 3502.72 (8000) |
|    |   | 32) CNA+semilag+dc: 3175.01 (8000) |
|    |   | 33) CNA+upwind+dc: 3330.177 (8000) |
|    |   | **Ra = 1e6** (Diff+Adv+Momentum)** |
|    |   | 34) Explicit+semilag+dc: 3327.062 (8000) |
|    |   | 35) CNA+semilag+dc: 3387.836 (8000) |
|    |   | 36) CNA+upwind+dc: 3277.684 (8000) |
| 13 | Blankenbach Benchmark | **Resolution: 50×50** |
|    |   | Ra = 1e4: 176.532 s (2108) |
|    |   | Ra = 1e5: 213.955 s (2406) |
|    |   | Ra = 1e6: 1003.409 s (8000) |
|    |   | **Resolution: 100×100** |
|    |   | Ra = 1e4: 3309.393 s (8000) |
|    |   | Ra = 1e5: 3232.410 s (5045) |
|    |   | Ra = 1e6: 5008.271 s (8000) |
|    |   | **Resolution Tests:** |
|    |   | Ra = 1e4: 5217.55 s |
|    |   | Ra = 1e5: 20459.1 s |
|    |   | Ra = 1e6: 554874 s |



| Exercise Number | Type                               | Total Runtime [s] |
| --------------- | ---------------------------------- | ----------------- |
| 1               | Euler Advection                    | 5.741             |
| 2               | 1-D Heat Diffusion <br> (explicit) | 1) 9.342 (1.14) \ 2) 9.975 (1.087) \newline 3) 9.616 (1.017) |
| 3               | 1-D Heat Diffusion \n (implicit) | 1) 12.641 (1.749) <br/> 2) 13.537 (1.720) <br/> 3) 12.157 (1.731) | 
| 4               | 2-D Steady State Heat Equation     | 13.471 (2.069) |
| 5a              | 2-D Time-dep. Heat Equation (Plume) | 1) Explicit: 18.632 (4.799) <br/> 2) Implicit: 78.201 (60.828) |
| 5b              | 2-D Time-dep. Heat Equation (Sill)  | 1) Explicit: 17.856 (4.566) <br/> 2) Implicit: 39.626 (22.679) |
| 6               | 1-D Advection Scheme                | 1) FTCS: 11.644 (4.019) <br/> 2) Upwind: 11.03 (4.042) <br/> 3) Downwind: 11.465 (3.343) <br/> 4) LAX: 11.272 (3.405) <br/> 5) SLF: 11.312 (3.48) <br/> 6) Semilag: 14.318 (3.604) <br/> 7) Tracers: 15.662 (4.944) |
| 7               | 2-D Energy Equation                 | 1) Upwind+explicit: 17.468 (2.979) <br/> 2) Upwind+implicit: 31.441 (10.842) <br/> 3) Upwind+CNA: 32.065 (12.222) <br/> 4) Upwind+ADI: 121.665 (99.807) <br/> 5) Upwind+dc: 30.45 (11.401) <br/> 6) SLF+explicit: 22.002 (2.872) <br/> 7) SLF+implicit: 28.751 (8.911) <br/> 8) SLF+CNA: 29.876 (9.28) <br/> 9) SLF+ADI: 120.018 (98.74) <br/> 10) SLF+dc: 28.179 (9.379) <br/> 11) SL+explicit: 22.969 (6.334) <br/> 12) SL+implicit: 30.634 (9.682) <br/> 13) SL+CNA: 32.415 (11.566) <br/> 14) SL+ADI: 133.871 (110.704) <br/> 15) SL+dc: 35.364 (12.301) <br/> 16) Tracers+explicit: 89.499 (65.913) <br/> 17) Tracers+implicit: 108.464 (84.209) <br/> 18) Tracers+CNA: 112.171 (89.219) <br/> 19) Tracers+ADI: 186.71 (161.611) <br/> 20) Tracers+dc: 115.737 (89.089) |
| 8               | 1-D Stokes Equation                 | 8.489 (0.212) |
| 9               | 2-D Falling Block (steady state)    | 20.764 (0.891) | 
| 10              | 2-D Falling Block (time-dep.)       | 1) Upwind: 33.295 (3.68) <br/> 2) SLF: 32.333 (3.629) <br/> 3) SL: 37.537 (3.755) <br/> 4) Tracers: 36.291 (6.12) |
| 11<sup>1</sup> | Thermal Convection (dim) <br/> **Resolution: 150x50**    | **Ra = 1e4** (Diff+Adv+Momentum) <br/> 1) Explicit+upwind+direct: 262.912 (1703) <br/> 2) Implicit+upwind+direct: 283.031 (1709) <br/> 3) CNA+upwind+direct: 284.944 (1708) <br/> 4) ADI+upwind+direct: 645.181 (1708) <br/> 5) DC+upwind+direct: 294.1 (1710) <br/> 6) Explicit+SLF+direct: 354.35 (2341) <br/> 7) Implicit+SLF+direct: 382.307 (2345) <br/> 8) CNA+SLF+direct: 386.764 (2344) <br/> 9) ADI+SLF+direct: 892.411 (2344) <br/> 10) DC+SLF+direct: 492.168 (2346) <br/> 11) Explicit+semilag+direct: 334.941 (1743) <br/> 12) Implicit+semilag+direct: 367.501 (1746) <br/> 13) CNA+semilag+direct: 371.289 (1744) <br/> 14) ADI+semilag+direct: 704.852 (1744) <br/> 15) DC+semilag+direct: 390.767 (1749) <br/> 16) Explicit+upwind+DC: 371.091 (1703) <br/> 17) Implicit+upwind+DC: 358.323 (1709) <br/> 18) CNA+upwind+DC: 338.862 (1708) <br/> 19) ADI+upwind+DC: 654.069 (1708) <br/> 20) DC+upwind+DC: 316.456 (1710) <br/> 21) Explicit+SLF+DC: 394.471 (2341) <br/> 22) Implicit+SLF+DC: 417.712 (2345) <br/> 23) CNA+SLF+DC: 412.717 (2344) <br/> 24) ADI+slf+DC: 890.714 (2344) <br/> 25) DC+SLF+DC: 403.217 (2346) <br/> 26) Explicit+semilag+DC: 286.046 (1743) <br/> 27) Implicit+semilag+DC: 295.961 (1746) <br/> 28) CNA+semilag+DC: 288.067 (1744) <br/> 29) ADI+semilag+DC: 671.159 (1744) <br/> 30) DC+semilag+DC: 357.894 (1749) <br/> <br/> **Ra = 1e5** (Diff+Adv+Momentum) <br/> 1) Explicit+semilag+dc: 489.869 (2433) <br/> 2) CNA+semilag+dc: 515.539 (2448) <br/> 3) CNA+upwind+dc: 1327.856 (8000) <br/><br/> **Ra = 1e6** (Diff+Adv+Momentum) <br/> 4) Explicit+semilag+dc: 1280.077 (8000) <br/> 5) CNA+semilag+dc: 1324.29 (8000) <br/> 6) CNA+upwind+dc: 1272.997 (8000)  |
| 12            | Thermal Convection (scaled) <br/> **Resolution: 150x50**    | **Ra = 1e4** (Diff+Adv+Momentum) <br/> 1) Explicit+upwind+direct: 310.208 (1950) <br/> 2) Implicit+upwind+direct: 333.82 (1957) <br/> 3) CNA+upwind+direct: 331.791 (1955) <br/> 4) ADI+upwind+direct: 739.821 (1955) <br/> 5) DC+upwind+direct: 352.163 (1953) <br/> 6) Explicit+SLF+direct: 437.435 (2704) <br/> 7) Implicit+SLF+direct: 468.642 (2708) <br/> 8) CNA+SLF+direct: 464.97 (2707) <br/> 9) ADI+SLF+direct: 1013.196 (2707) <br/> 10) DC+SLF+direct: 476.063 (2706) <br/> 11) Explicit+semilag+direct: 358.84 (1999) <br/> 12) Implicit+semilag+direct: 341.811 (2002) <br/> 13) CNA+semilag+direct: 384.29 (2000) <br/> 14) ADI+semilag+direct: 754.453 (2000) <br/> 15) DC+semilag+direct: 359.703 (2001) <br/> 16) Explicit+upwind+DC: 683.465 (1950) <br/> 17) Implicit+upwind+DC: 614.899 (1957) <br/> 18) CNA+upwind+DC:  649.386 (1955) <br/> 19) ADI+upwind+DC: 1048.092 (1955) <br/> 20) DC+upwind+DC:  640.199 (1953) <br/> 21) Explicit+SLF+DC:  835.568 (2704) <br/> 22) Implicit+SLF+DC:  862.991 (2708) <br/> 23) CNA+SLF+DC:  868.07 (2707) <br/> 24) ADI+slf+DC:  1402.997 (2707) <br/> 25) DC+SLF+DC: 876.505 (2706) <br/> 26) Explicit+semilag+DC: 613.413 (1999) <br/> 27) Implicit+semilag+DC: 652.126 (2002) <br/> 28) CNA+semilag+DC: 638.264 (2000) <br/> 29) ADI+semilag+DC:  1027.153 (2000) <br/> 30) DC+semilag+DC: 620. 945 (2001) <br/> <br/> **Ra = 1e5** (Diff+Adv+Momentum) <br/> 31) Explicit+semilag+dc:  3502.72 (8000) <br/>  32) CNA+semilag+dc: 3175.01 (8000) <br/> 33) CNA+upwind+dc: 3330.177 (8000) <br/><br/> **Ra = 1e6** (Diff+Adv+Momentum) <br/> 34) Explicit+semilag+dc: 3327.062 (8000) <br/> 35) CNA+semilag+dc: 3387.836 (8000) <br/> 36) CNA+upwind+dc: 3277.684 (8000) <br/> | 
| 13            | Blankenbach Benchmark         | **Resolution: 50x50** <br/><br/> 1) Ra = 1e4: 176.532 s(2108) <br/> 2) Ra = 1e5: 213.955 s (2406) <br/> 3) Ra = 1e6: 1003.409 s (8000) <br/> <br/> **Resolution: 100x100** <br/><br/> 1) Ra = 1e4: 3309.393 s (8000) <br/> 2) Ra = 1e5: 3232.410 s (5045) <br/> 3) Ra = 1e6: 5008.271 s (8000) <br/><br/> **Resolution Tests:** <br/><br/> Ra = 1e4: 5217.55 s <br/> Ra = 1e5: 20459.1 s <br/> Ra = 1e6: 554874 s |

---

¹: For the thermal convection, only the runtime where compilation and allocation overhead is avoided are shown, for the sake of simplicity. The number of iterations is shown in the parentheses. For the higher $Ra$ number models only the fastest and most accurate (see [here](../man/examples/GaussianDiffusion2D.md)) combinations of solvers from the low $Ra$ case haven been run. 



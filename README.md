# GeoModBox.jl
**Geodynamic Modelling** Tool**Box** is a julia package mainly used for teaching purposes. The package provides different finite differences, fully staggered, discretization schemes for the governing equations to describes a geodynamic problem. The geoverning equations are the conservation equations of **energy**, **mass**, **momentum**, and **compositon**. 

<!-- I will add some question directly here for the time beeing: I would like to add similar informations to this package as I did in my [FDCSGm repository](https://github.com/LukasFuchs/FDCSGm). However, I am not familiar with the doc option in github yet, that's why I added all the information in a README.md file for each module etc. I guess the doc option would be more suitable to add the details and simply add some figures and general information in the README.md. Have you ever worked with doc in github? -->

------------------
------------------
## Staggered Finite Difference 
- brief general description of a staggered finite difference scheme
- Energy equation
- Momentum equation
------------------
------------------
## Energy equation 
&emsp; The conservation of energy is a fundamental principle in physics and defines that the loss and generation of energy needs to be equal. In terms of a geodynamical problem, energy can be described by temperature, which is transported mainly through *conductive* and *convective* processes, such that a general energy equation is defined as followed (assuming only radioactive heat sources):

$$
(\frac{\partial E}{\partial t} + \overrightarrow{v} \cdot \nabla E) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H, \tag{1}
$$

where the energy is described as $E=c_{p} \rho T$, and *c<sub>p</sub>* is the specific heat capacity [J/kg/K], *ρ* is a reference density [kg/m<sup>3</sup>], *T* is the temperature [K], *t* is the time [s], $\overrightarrow{v}$ is the velocity vector [m/s], *q<sub>i</sub>* is the heat flux in direction of *i*  [W/m<sup>2</sup>], *∂/∂xi* is a directional derivative in direction of *i*, and *H* the heat production rate per mass [W/kg]. The repeated index means a summation of derivatives. This conservation law contains the variation of the heat flux in a certain direction, where the heat flux is defined by the Fourier’s law as followed: 

$$
\overrightarrow{q} = - k \nabla T, \tag{} 
$$

where *k* is the thermal conductivity [W/m/K]. The heat flux is the amount of heat that passes through a unit surface area, per unit time and is positive in the direction of decreasing temperature, that is in the case when the temperature gradient is negative. The *temperature conservation equation* in an Eulerian form can then be written as: 

$$
\rho c_p (\frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T) = -\frac{\partial q_i}{\partial x_i} + \rho H. \tag{}
$$

- General formulation of the energy equation (advection + diffusion term + radiogenic heating)
- Operator splitting
- Solution of the conductive part
------------------
### Heat Diffusion
- Conducitve part of the energy equation + radiogenic heating + ?
- Constant thermal parameters
- Variable thermal paramters
------------------
#### Numerical Schemes

All numerical schemes methods can be used in the [thermal convection code]() and the [Blankenbach Benchmark]() and are generally available to chose in the code. A more detailed analysis on the accuracy of each discretization scheme and the effect of the grid resolution is given in the [Gaussian Diffusion Benchmark](). 

##### Explicit, FTCS
##### Implicit, FTCS
##### Cranck-Nicolson Approach (CNA)
##### Alternating Direction Implicit (ADI)
------------------
##### Thermal Boundary Conditions
1. Dirichlet
2. Neumann
------------------
#### Steady State Solution
##### Poisson solution, constant *k*
##### Poisson solution, variable *k*
------------------
------------------
## Advection Equation
### Discretization Schemes
#### Upwind 
#### Staggered Leap Frog (SLF)
#### Semi-Lagrangian
#### Passive Tracers
------------------
------------------
## Stokes Equation
------------------
------------------
## Benchmarks
------------------
------------------
## Examples
------------------
------------------

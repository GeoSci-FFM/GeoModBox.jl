# Momentum Conservation Equation

## General Information

&emsp;On geological time scales, Earth’s mantle and lithosphere do behave like a fluid and move and deform. A fluid does generally move due to forces acting on it whereas the forces must be in balance. In general, there are three major forces one might consider, i.e., inertia, surface, and volumetric forces. A common equation to describes such motion is given by: 

$$ 
\begin{equation}
\rho \frac{D \overrightharpoon{v}}{Dt} = \nabla \cdot \boldsymbol{\sigma} + \rho \boldsymbol{g},
\end{equation}
$$

where $\rho$ is the density [kg/m<sup>3</sup>], $\overrightharpoon{v}$ is the velocity vector [m/s], $\boldsymbol{\sigma}$ is the *Cauchy stress tensor* [Pa], $\boldsymbol{g}$ is the gravitational acceleration [m/s<sup>2</sup>], and the term on the left-hand side is the Lagrangian time derivative which is in Eulerian form $\frac{D}{Dt} = \frac{\partial{}}{\partial{t}} + \overrightharpoon{v} \cdot \nabla$. 

The *Cauchy stress tensor* is given by: 

$$
\begin{equation}
\boldsymbol{\sigma} = -\nabla{P} + \nabla \cdot \boldsymbol{\tau},
\end{equation}
$$

where *P* is the total pressure (*P = P<sub>dynamic</sub> + P<sub>hydrostatic</sub>*) and $\boldsymbol{\tau}$ the *deviatoric stress tensor*. 

In Eulerian form, equation $(1)$ is given by (**Navier-Stokes equation**):

$$
\begin{equation}
\rho \left(\frac{\partial{v_{i}}}{\partial{t}} + v_{j}\frac{v_{i}}{\partial{x_{j}}}\right) = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\tau_{ij}}{\partial{x_j}} + \rho g_{i},
\end{equation}
$$

where summation over repeated indices is implied.

### Constitutive Relation

&emsp;To solve equation $(3)$, one needs to define a rheology which, for a purely viscous medium, can be given by a constitutive relationship between stress and strain rate in the form of, e.g.:

$$
\begin{equation}
\tau_{ij} = 2 \eta \cdot \dot{\varepsilon}_{ij},
\end{equation}
$$

where $\eta$ is the dynamic viscosity in [Pa s] and $\dot{\varepsilon}_{ij}$ the *strain rate tensor* in [1/s] and given by: 

$$
\begin{equation}
\dot{\varepsilon}_{ij} = \frac{1}{2} \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right),
\end{equation}
$$

### Stokes Equation

&emsp;Assuming that the inertia forces are negligible in comparison to the gravitational forces, one can further simplify equation $(3)$ to:

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\tau_{ij}}{\partial{x_j}} + \rho g_{i},
\end{equation}
$$

or in the form of the unknowns $v_x$, $v_y$, and $P$:

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial}{\partial{x_j}} \eta \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right) + \rho g_{i}.
\end{equation}
$$

Assuming a constant viscosity further simplifies equation $(7)$ to: 

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \eta \left(\frac{\partial^2{v_i}}{\partial{x_j^2}} + \frac{\partial^2{v_j}}{\partial{x_i^2}}\right) + \rho g_{i}.
\end{equation}
$$

### Continuum Equation

&emsp;Equation $(8)$ provides us two equations for our three unknowns. Thus, one needs to also consider the mass conservation equation (we do work with a continuum), where one can further simplify the problem by assuming an incompressible fluid (i.e., Boussinesq-approximation):

$$
\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}
$$

Equations $(8)$ and $(9)$ enable us to solve for the three unknowns $v_x$, $v_y$, and $P$. 

### Equation of State

&emsp;The buoyance term on the right-hand side of equation $(7)$, that is the density term which is temperature dependent (and pressure, but I do neglect this effect here so far), can be approximated with the so-called *equation of state* for the density. Here, its is a linear approximation of the change of density due to temperature variations and can be defined as:

$$
\begin{equation}
\rho = \rho_0 (1-\alpha T),
\end{equation}
$$

where *ρ<sub>0</sub>* is the reference density and *α* the thermal expansion coefficient [1/K]. 

<!--
- Scaling 
 -->

------------------
------------------

## Examples
<!--
### Channel Flow
&emsp; Assuming the horizontal pressure gradient is constant and flow within a channel is only driven by the pressure and/or by a constant horizonal velocity at the surface (or at the bottom, or both), the stokes equation describes the horizontal flow velocity within the channel and simplifies to: 

$$
\frac{\partial P}{\partial x} = \frac{\partial \tau_{x,z}}{\partial z}
$$,

where *P* is the pressure and $\tau_{x,z}$ is the deviatoric shear stress, which is defined as: 

$\tau_{x,z} = 2 \eta \dot{\varepsilon}_{x,z} = \eta \frac{\partial v_x}{\partial z}$,

where $\eta$ is the dynamic viscosity and $\dot{\varepsilon}_{x,z}$ is the deviatoric shear strain-rate, which is defined as: 

$\dot{\varepsilon}_{x,z} = \frac{1}{2}(\frac{\partial v_x}{\partial z} + \frac{\partial v_z}{\partial x})$.

For the given setup I can assume that the vertical velocity is zero and thus equation (3) simplifies to the last expression of equation (2).

&emsp; This directory contains a script to calculate the horizontal velocity for a two-dimensional Couette(-Poiseuille) channel flow with constant and logarithmically, with depth varying viscosity and to compare the numerical solution with its analytical solution. The depth-dependent viscosity is defined as: 

$\eta = \eta_0 exp(log(m) \frac{H-z}{H})$,

where *m* is the viscosity ratio of $\frac{\eta_1}{\eta_0}$, $\eta_0$ and $\eta_1$ are the bottom and surface viscosities, respectively, *H* is the model height, and *z* the depth. 

&emsp;Considering the definition of the viscosity as given in equation (4), one can derive an analytical solution of the horizontal velocity from the 1-D stokes equation in *x*-direction by twice integrating equation (1). The analytical solution with depth depends on the viscosity ratio *m*, the horizontal pressure gradient $\frac{\partial P}{\partial x}$, and the shear velocity at the surface $v_{x,0}$. For an upward pointing coordinate system (*z* positive) the analytical solution is given as: 

$v_{x,ana} =-\frac{1}{2 \eta_0} \frac{\partial P}{x} (Hz - z^2) + v_{x,0}\frac{z}{H}$,&emsp;&emsp; if $m = 1$, and &emsp;&emsp;&emsp; (5)

$v_{x,ana} = -\frac{\partial P}{\partial x} \frac{H}{\eta_0 log(m)} (\frac{m^{-\frac{z}{H}}}{m-1}(z(m-1)+H) - \frac{H}{m-1})-m^{-\frac{z}{H}} m \frac{v_{x,0}}{m-1} + \frac{v_{x,0}m}{m-1}$, &emsp;&emsp; if $m \neq 0$.&emsp;&emsp;&emsp; (6)

&emsp;The numerical solution is calculated using fixed boundary velocities, which are defined by the analytical solution of the horizontal velocity as defined in equations (5) and (6) and I simply flip the analytical solution so that it fits to the downward point coordinate system I use in the code. -->

---------------------------------------------------------------

<!-- Channel flow example in quasi 2-D:
    - prescribe velocity boundary conditions by the analytical solution of a channel flow 
-->
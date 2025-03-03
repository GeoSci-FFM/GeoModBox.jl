# Energy Conservation Equation

&emsp; The conservation of energy is a fundamental principle in physics and defines that the loss and generation of energy needs to be equal within a closed system. In terms of a geodynamical problem, energy is described by temperature, which is transported mainly through *conductive* and *convective* processes, such that a general energy equation is defined as followed (assuming only radioactive heat sources):

$$
\begin{equation}
\left(\frac{\partial E}{\partial t} + \overrightarrow{v} \cdot \nabla E\right) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H,
\end{equation}
$$

where the energy is described as $E=c_{p} \rho T$, and *c<sub>p</sub>* is the specific heat capacity [J/kg/K], *ρ* is a reference density [kg/m<sup>3</sup>], *T* is the temperature [K], *t* is the time [s], $\overrightarrow{v}$ is the velocity vector [m/s], *q<sub>i</sub>* is the heat flux in direction of *i*  [W/m<sup>2</sup>], *∂/∂xi* is a directional derivative in direction of *i*, and *H* the heat production rate per mass [W/kg]. The repeated index means a summation of derivatives. This conservation law contains the variation of the heat flux in a certain direction, where the heat flux is defined by the Fourier’s law as followed: 

$$
\begin{equation}
\overrightarrow{q} = - k \nabla T,
\end{equation}
$$

where *k* is the thermal conductivity [W/m/K]. The heat flux is the amount of heat that passes through a unit surface area, per unit time and is positive in the direction of decreasing temperature, that is in the case when the temperature gradient is negative. The *temperature conservation equation* in an Eulerian form can then be written as: 

$$
\begin{equation}
\rho c_p \left(\frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla T\right) = -\frac{\partial q_i}{\partial x_i} + \rho H.
\end{equation}
$$

&emsp;This form of the temperature equation describes the variation of temperature due to a *conductive* (right hand side of the equation) and *convective* (left hand side of the equation) process. For a matter of simplicity, one can consider those terms in a separate manner and solve the *temperature conservation equation* using an *operator splitting* method, that is one first solves the *convective* part, followed by the *conductive* part. 

&emsp; This directory focusses on examples for the **conductive** part of the *temperature conservation equation* using [different numerical finite difference schemes](../../src/HeatEquation/) applicable to [one-](./1D/) and [two-](./2D/)dimensional problem. 

<!--
- Scaling
 -->
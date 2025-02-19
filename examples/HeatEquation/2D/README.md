# General Information 

## Energy Equation 

## Discretization 

### Explicit, FTCS (or Forward Euler Method)

$$
\begin{equation}
\frac{T_{i,j}^{n+1} - T_{i,j}^{n} }{\Delta t} = \kappa \left( \frac{T_{i-1,j}^{n} - 2T_{i,j}^{n} + T_{i+1,j}^{n}}{\Delta{x}^2} + \frac{T_{i,j-1}^{n} - 2T_{i,j}^{n} + T_{i,j+1}^{n}}{\Delta{z^2}} \right) + \frac{Q_{i,j}^n}{\rho c_p}
\end{equation}
$$

$$
\begin{equation}
T_{i,j}^{n+1} = T_{i,j}^{n} + s_x\left(T_{i-1,j}^{n} - 2T_{i,j}^{n} + T_{i+1,j}^{n}\right) + s_z \left(T_{i,j-1}^{n} - 2T_{i,j}^{n} + T_{i,j+1}^{n}\right) + \frac{Q_{i,j}^n \Delta{t}}{\rho c_p}
\end{equation}
$$

#### Boundary Conditions

## Numerical Schemes

### Implicit, FTCS (or Backward Euler Method)

#### Boundary Conditions

### Defection Correction Method

### Cranck-Nicolson Approach (CNA)

---------------
---------------

## Examples 

### Gaussian Diffusion ([]())

### Resolution Test Poisson Problem
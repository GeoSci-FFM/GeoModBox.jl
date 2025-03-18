# General Information 

<!-- What is this directory about-->

--------------------
--------------------

## Stokes Equation (1D)

*x-component*

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x}} + \frac{\partial{\tau_{xy}}}{\partial{y}}
\end{equation}
$$

## Discretization 

<!-- add figure of grid -->

***Constant Viscosity***

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x}} + \eta\frac{\partial^2{v_x}}{\partial{y^2}}.
\end{equation}
$$

<!--Approximation of partial derivatives with finite difference operators ... -->

$$
\begin{equation}
\frac{\partial{P}}{\partial{x}}=av_{x,j-1}+bv_{x,j}+cv_{x,j+1}, 
\end{equation}
$$

where

$$
a = c = \frac{\eta}{\Delta{y^2}},\ \textrm{and}\ b = -\frac{2\eta}{\Delta{y^2}}.
$$

***Variable Viscosity***

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x}} + \frac{\partial}{\partial{y}}\left(\eta\frac{\partial{v_x}}{\partial{y}}\right)
\end{equation}
$$

$$
\begin{equation}
\frac{\partial{P}}{\partial{x}}=\frac{\eta_{j+1}\frac{\partial{v_x}}{\partial{y}}\vert_{j+1}-\eta_{j}\frac{\partial{v_x}}{\partial{y}}\vert_{j}}{\Delta{y}}
\end{equation}
$$

$$
\begin{equation}
\frac{\partial{P}}{\partial{x}}=\frac{\eta_{j+1}\frac{v_{x,j+1}-v_{x,j}}{\Delta{y}}-\eta_{j}\frac{v_{x,j}-v_{x,j-1}}{\Delta{y}}\vert_{j}}{\Delta{y}}
\end{equation}
$$

$$
\begin{equation}
\frac{\partial{P}}{\partial{x}}=av_{x,j-1}+bv_{x,j}+cv_{x,j+1}, 
\end{equation}
$$

where

$$
\begin{equation}
a = \frac{\eta_j}{\Delta{y^2}}, b = -\frac{\eta_j+\eta_{j+1}}{\Delta{y^2}},\ \textrm{and}\ c = \frac{\eta_{j+1}}{\Delta{y^2}}
\end{equation}
$$

### Boundary Conditions

**Dirichlet**

*Bottom* (j=1)

$$\begin{equation}
V_{G,S} = 2V_{BC,S} - v_{x,1}
\end{equation}$$

*Top* (j=nc)

$$\begin{equation}
V_{G,N} = 2V_{BC,N} - v_{x,nc}
\end{equation}$$

**Neumann**

*Bottom* (j=1)

$$\begin{equation}
V_{G,S} = v_{x,1} - c_s\Delta{y},
\end{equation}$$

*Top* (j=nc)

$$\begin{equation}
V_{G,N}=v_{x,nc} + c_N\Delta{y},
\end{equation}$$

where 

$$\begin{equation}
c_S = \frac{dv_x}{dy}\vert_{j=1},\ \textrm{and}\ c_N=\frac{dv_x}{dy}\vert_{j=nc}
\end{equation}$$

<!-- Change of the coefficient matrix and the right hand side ... -->

**Dirichlet**

*Bottom* (j=1)

$$\begin{equation}
\left(b-a\right)v_{x,1}+cv_{x,2} = \frac{\partial{P}}{\partial{x}} - 2aV_{BC,S}
\end{equation}
$$

*Top* (j=nc)

$$\begin{equation}
av_{x,nc-1}+\left(b-c\right)v_{x,nx} = \frac{\partial{P}}{\partial{x}} - 2cV_{BC,N}
\end{equation}
$$

**Neumann**

*Bottom* (j=1)

$$\begin{equation}
\left(b+a\right)v_{x,1}+cv_{x,2} = \frac{\partial{P}}{\partial{x}} + ac_SΔy
\end{equation}
$$

*Top* (j=nc)

$$\begin{equation}
av_{x,nc-1}+\left(b+c\right)v_{x,nx} = \frac{∂P}{∂x} - ac_NΔy
\end{equation}
$$

### Solution 

#### Direct

<!-- rhs includes the boundary information! -->

$$\begin{equation}
\bold{K} \cdot v_x = rhs
\end{equation}
$$

$$\begin{equation}
v_x = \bold{K} ∖ rhs
\end{equation}
$$

#### Defect Correction 

$$\begin{equation}
R = -\frac{∂P}{∂x} + \frac{∂τ_{xy}}{∂y}
\end{equation}
$$

$$\begin{equation}
R = -\frac{∂P}{∂x} + \bold{K} \cdot v_x
\end{equation}
$$

$$\begin{equation}
R_i = -\frac{∂P}{∂x} + \bold{K_i} \cdot v_{x,i}
\end{equation}
$$s

$$\begin{equation}
0 = -\frac{∂P}{∂x} + \bold{K}\left(v_{x,i}+ \delta{v_x} \right) = \bold{K}\cdot v_{x,i} -\frac{∂P}{∂x} + \bold{K}\cdot \delta{v_x} = R_i + \bold{K} \cdot \delta{v_x}
\end{equation}
$$

$$\begin{equation}
R_i = -\bold{K}\cdot{\delta{v_x}}
\end{equation}
$$

$$\begin{equation}
\delta{v_x} = -\bold{K}^{-1}R_i
\end{equation}
$$

$$\begin{equation}
v_x^n = v_{x,i} + \delta{v_x}
\end{equation}
$$

--------------------
--------------------

## Examples

- [Channel Flow (using the defect correction method)](./ChannelFlow_1D.jl)
<!-- 
- 1D case 
-- Discretized equations
-- Solving the equations 
--- Direct solution 
--- Defection corrections solution
- 2D case 
-->
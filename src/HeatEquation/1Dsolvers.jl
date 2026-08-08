using ExtendableSparse
# ======================================================================= #
# Time-dependent solvers, constant thermal parameters =================== #
# ======================================================================= #
"""
    ForwardEuler1Dc!(explicit, κ, Δx, Δt, nc, BC;
                     Q=zeros(nc), ρ=3200.0, cp=1200.0)

Solves the one-dimensional transient heat equation using an explicit
Forward Euler finite-difference scheme.

Temperature is defined at cell centroids, while the heat flux is evaluated
at the cell boundaries. Ghost nodes are used to impose Dirichlet or Neumann
boundary conditions. Thermal diffusivity is assumed to be constant.

Optional internal heating can be included through the volumetric heat
production term `Q`.

# Arguments

    explicit    : Structure or tuple containing:
                  `T`, the temperature array on the centroids, and
                  `T_ex`, the extended temperature array including ghost nodes.
    κ           : Thermal diffusivity.
    Δx          : Grid spacing.
    Δt          : Time step.
    nc          : Number of centroid nodes.
    BC          : Structure or tuple defining the boundary-condition types
                  and values at the western and eastern boundaries.

# Keyword Arguments

    Q           : Volumetric heat production rate defined on the centroids
                  (default: `zeros(nc)`).
    ρ           : Density (default: `3200.0`).
    cp          : Specific heat capacity (default: `1200.0`).

# Notes

The temperature update is explicit and second-order accurate in space but
first-order accurate in time. The timestep must therefore satisfy the
diffusive stability condition.

After each timestep, the updated centroid temperatures are copied back to
the interior of `explicit.T_ex`. The ghost-node values are updated again 
when the routine is called for the next timestep.
"""
function ForwardEuler1Dc!( explicit, κ, Δx, Δt, nc, BC; 
                                Q = zeros(nc), ρ=3200.0, cp=1200.0 )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    # Define boundary conditions ---------------------------------------- #
    # West ---
    explicit.T_ex[1]    =   (BC.type.W==:Dirichlet) * (2 * BC.val.W - explicit.T_ex[2]) + 
                            (BC.type.W==:Neumann) * (explicit.T_ex[2] - BC.val.W*Δx)
    # East --
    explicit.T_ex[end]  =   (BC.type.E==:Dirichlet) * (2 * BC.val.E - explicit.T_ex[nc+1]) +
                            (BC.type.E==:Neumann) * (explicit.T_ex[nc+1] + BC.val.E*Δx)
    for i = 1:nc
        # Calculate temperature at point i for the new time ---        
        explicit.T[i] =   explicit.T_ex[i+1] + κ * Δt * 
                (explicit.T_ex[i + 2] - 2.0 * explicit.T_ex[i+1] + explicit.T_ex[i]) / Δx^2 +
                Q[i] * Δt / ρ / cp
    end
    explicit.T_ex[2:end-1]  .=  explicit.T
end

"""
    ComputeResiduals1Dc!(
        R, T, T_ex, T0, T_ex0, ∂2T, κ, BC, Δx, Δt;
        C=0.0, Q=0.0, ρ=3200.0, cp=1200.0
    )

Computes the residual of the one-dimensional transient heat equation for
constant thermal parameters.

Temperature is defined at cell centroids, while second spatial derivatives
are evaluated using central finite differences and ghost nodes. Dirichlet
and Neumann boundary conditions are supported at the western and eastern
boundaries.

The temporal discretization is controlled by `C` and can represent Backward
Euler, Crank–Nicolson, or Forward Euler time integration. Optional internal
heating can be included through the volumetric heat-production term `Q`.

# Arguments

    R        : Residual vector defined on the centroids.
    T        : Temperature at the current iteration or new time level.
    T_ex     : Extended current-temperature array including ghost nodes.
    T0       : Temperature at the previous time level.
    T_ex0    : Extended previous-temperature array including ghost nodes.
    ∂2T      : Structure or tuple containing the second spatial derivatives
               `∂x2` and `∂x20`.
    κ        : Thermal diffusivity.
    BC       : Structure or tuple defining the boundary-condition types and
               values at the western and eastern boundaries.
    Δx       : Grid spacing.
    Δt       : Time step.

# Keyword Arguments

    C        : Temporal weighting parameter (default: `0.0`):
               `0.0` for Backward Euler,
               `0.5` for Crank–Nicolson,
               `1.0` for Forward Euler.
    Q        : Volumetric heat-production rate (default: `0.0`).
    ρ        : Density (default: `3200.0`).
    cp       : Specific heat capacity (default: `1200.0`).

# Notes

The residual is evaluated as a weighted combination of the diffusion term at
the previous and current time levels. Backward Euler and Crank–Nicolson require
the current-temperature contribution, whereas Forward Euler uses only the
previous temperature field.

Backward Euler is first-order accurate in time, Crank–Nicolson is second-order
accurate in time, and Forward Euler is first-order accurate in time. The
spatial discretization is second-order accurate.
"""
function ComputeResiduals1Dc!( R, T, T_ex, T0, T_ex0, ∂2T, κ, BC, Δx, Δt;
                C=0,Q=0.0,ρ=3200.0,cp=1200.0)
    if C < 1
        # Implicit
        T_ex[2:end-1]   .=  T    
        T_ex[1]         =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex[2]) + (BC.type.W==:Neumann)*(T_ex[2] - BC.val.W*Δx)
        T_ex[end]       =   (BC.type.E==:Dirichlet)*(2*BC.val.E - T_ex[end-1]) + (BC.type.E==:Neumann)*(T_ex[end-1] + BC.val.E*Δx)
        @. ∂2T.∂x2      =   (T_ex[3:end] - 2 * T_ex[2:end-1] + T_ex[1:end-2]) / Δx^2
        if C==0.5
            # CNA
            T_ex0[2:end-1]  .=  T0    
            T_ex0[1]        =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex0[2]) + (BC.type.W==:Neumann)*(T_ex0[2] - BC.val.W*Δx)
            T_ex0[end]      =   (BC.type.E==:Dirichlet)*(2*BC.val.E - T_ex0[end-1]) + (BC.type.E==:Neumann)*(T_ex0[end-1] + BC.val.E*Δx)
            @. ∂2T.∂x20     =   (T_ex0[3:end] - 2 * T_ex0[2:end-1] + T_ex0[1:end-2]) / Δx^2
        end
    else
        # explicit
        T_ex0[2:end-1]  .=   T0    
        T_ex0[1]        =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex0[2]) + 
                                (BC.type.W==:Neumann)*(T_ex0[2] - BC.val.W*Δx)
        T_ex0[end]      =   (BC.type.E==:Dirichlet)*(2*BC.val.E - T_ex0[end-1]) + 
                                (BC.type.E==:Neumann)*(T_ex0[end-1] + BC.val.E*Δx)
        @. ∂2T.∂x20     =   (T_ex0[3:end] - 2 * T_ex0[2:end-1] + T_ex0[1:end-2]) / Δx^2
    end
    # Calculate residual ------------------------------------------------ #
    @. R             =   (T - T0)/Δt - κ*((1-C)*∂2T.∂x2 + C*∂2T.∂x20) - Q/ρ/cp
    # ------------------------------------------------------------------- #
end

"""
    AssembleMatrix1Dc!(κ, Δx, Δt, nc, BC, K; C=0.0)

Assembles the coefficient matrix for the one-dimensional transient heat
equation with constant thermal diffusivity.

Temperature is defined at cell centroids, and the spatial diffusion term is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the matrix
coefficients at the western and eastern boundaries.

The temporal discretization is controlled by `C`. For `C < 1`, the matrix
contains the implicit contribution of the diffusion term. For `C = 1`, the
diffusion contribution vanishes and the matrix reduces to the temporal term.

# Arguments

    κ        : Thermal diffusivity.
    Δx       : Grid spacing.
    Δt       : Time step.
    nc       : Number of centroid nodes.
    BC       : Structure or tuple defining the boundary-condition types and
               values at the western and eastern boundaries.
    K        : Coefficient matrix to be assembled in place.

# Keyword Arguments

    C        : Temporal weighting parameter (default: `0.0`):
               `0.0` for Backward Euler,
               `0.5` for Crank–Nicolson,
               `1.0` for Forward Euler.

# Notes

The matrix corresponds to the current-time contribution of the generalized
time-discretization scheme used by `ComputeResiduals1Dc!`. Backward Euler and
Crank–Nicolson require a diffusion contribution in the matrix, whereas Forward
Euler is fully explicit.

The assembled matrix is tridiagonal. Boundary conditions modify the diagonal
coefficients of the first and last equations. After assembly, `flush!(K)` is
called to finalize the sparse matrix.
    
"""
function AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K;C=0 )
    # Define coefficients ---
    a   =   (κ*(1-C)) / Δx^2
    b   =   1 / Δt
    # Loop over the grid points ---    
    for i = 1:nc  
        # Equation number ---
        ii          =   i
        # Stencil ---
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==nc ? false  : true
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        if inE
            K[ii,iE]    = - a
        end
        K[ii,iC]        =   (2 + DirW + DirE - NeuW - NeuE)*a + b
        if inW 
            K[ii,iW]    = - a
        end        
    end
    flush!(K)
end



"""
    BackwardEuler1Dc!(
        implicit, κ, Δx, Δt, nc, BC, K, rhs;
        Q=0.0, ρ=3200.0, cp=1200.0
    )

Solves the one-dimensional transient heat equation using an implicit Backward
Euler finite-difference scheme with constant thermal diffusivity.

Temperature is defined at cell centroids, while the diffusion term is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the coefficient
matrix and right-hand-side vector.

Optional internal heating can be included through the volumetric
heat-production term `Q`.

# Arguments

    implicit    : Structure or tuple containing `T`, the temperature array
                  defined on the centroids. The array is updated in place with
                  the solution at the new time level.
    κ           : Thermal diffusivity.
    Δx          : Grid spacing.
    Δt          : Time step.
    nc          : Number of centroid nodes.
    BC          : Structure or tuple defining the boundary-condition types and
                  values at the western and eastern boundaries.
    K           : Coefficient matrix for the linear system of equations.
    rhs         : Right-hand-side vector.

# Keyword Arguments

    Q           : Volumetric heat-production rate (default: `0.0`).
    ρ           : Density (default: `3200.0`).
    cp          : Specific heat capacity (default: `1200.0`).

# Notes

Backward Euler is first-order accurate in time and unconditionally stable for
the linear diffusion equation. The spatial discretization is second-order
accurate.

The coefficient matrix and right-hand-side vector are assembled inside the
routine. The temperature at the new time level is obtained by solving the
resulting linear system, and `implicit.T` is updated in place.
"""
function BackwardEuler1Dc!( implicit, κ, Δx, Δt, nc, BC , K, rhs; 
                                Q = 0.0, ρ = 3200.0, cp = 1200.0 )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    # Define coefficients ---
    a   =   κ / Δx^2
    b   =   1 / Δt
    # Multiply rhs with 1/Δt ---    
    @. rhs  =   b * implicit.T + Q/ρ/cp 
    # Loop over the grid points ---
    for i = 1:nc  
        # Equation number ---
        ii          =   i
        # Stencil ---
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==nc ? false  : true
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        if inE K[ii,iE]    = - a end
        K[ii,iC]        =   (2 + DirW + DirE - NeuW - NeuE)*a + b
        if inW K[ii,iW]    = - a end            
        # Modify right hand side due to boundary conditions ------------- #        
        # West boundary ---
        rhs[i]  +=  2*a*BC.val.W * DirW - 
                        a*BC.val.W*Δx * NeuW +
                        2*a*BC.val.E * DirE + 
                        a*BC.val.E*Δx * NeuE
    end
    # ------------------------------------------------------------------- #    
    # Calculate temperature at new time step ---------------------------- #
    implicit.T  .=   K \ rhs
    # ------------------------------------------------------------------- #    
end

"""
    CNA1Dc!(
        cna, κ, Δx, Δt, nc, BC, K1, K2;
        Q=0.0, ρ=3200.0, cp=1200.0
    )

Solves the one-dimensional transient heat equation using the Crank–Nicolson
finite-difference scheme with constant thermal diffusivity.

Temperature is defined at cell centroids, while the diffusion term is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the coefficient
matrices and right-hand-side vector.

The Crank–Nicolson method combines the diffusion operators evaluated at the
current and new time levels with equal weighting, resulting in a second-order
accurate temporal discretization. Optional internal heating can be included
through the volumetric heat-production term `Q`.

# Arguments

    cna         : Structure or tuple containing `T`, the temperature array
                  defined on the centroids. The array is updated in place with
                  the solution at the new time level.
    κ           : Thermal diffusivity.
    Δx          : Grid spacing.
    Δt          : Time step.
    nc          : Number of centroid nodes.
    BC          : Structure or tuple defining the boundary-condition types and
                  values at the western and eastern boundaries.
    K1          : Coefficient matrix associated with the unknown temperature
                  field at the new time level.
    K2          : Coefficient matrix associated with the known temperature
                  field at the previous time level.

# Keyword Arguments

    Q           : Volumetric heat-production rate (default: `0.0`).
    ρ           : Density (default: `3200.0`).
    cp          : Specific heat capacity (default: `1200.0`).

# Notes

The Crank–Nicolson scheme is second-order accurate in both space and time.
For linear diffusion problems it is unconditionally stable, although large
time steps may produce weak temporal oscillations.

The coefficient matrices `K1` and `K2` are assembled inside the routine. The
right-hand-side vector is constructed from the previous temperature field and
the volumetric heat-production term before solving the resulting linear
system. The computed temperature field overwrites `cna.T`.
"""
function CNA1Dc!( cna, κ, Δx, Δt, nc, BC, K1, K2; 
                        Q=0.0,ρ=3200.0,cp=1200.0 )
# ======================================================================= #
# LF; 19.09.2024 - Version 1.0 - Julia                                    #
# ======================================================================= #    
    rhs     = zeros(length(cna.T))
    # Define coefficients ---
    a       =   κ / 2 / Δx^2
    b       =   1 / Δt
    # Loop over the grid points ---
    for i = 1:nc
        # Equation number ---
        ii          =   i
        # Stencil ---
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==nc ? false  : true
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        if inE
            K1[ii,iE]   =   - a
            K2[ii,iE]   =   a            
        end
        K1[ii,iC]       =   b + (2 + DirW + DirE - NeuW - NeuE) * a
        K2[ii,iC]       =   b - (2 + DirW + DirE - NeuW - NeuE) * a 
        if inW 
            K1[ii,iW]   =   - a
            K2[ii,iW]   =   a
        end                    
    end
    # ------------------------------------------------------------------- #
    # Berechnung der rechten Seite -------------------------------------- #
    rhs     .=   K2 * cna.T .+ Q/ρ/cp
    # ------------------------------------------------------------------- #        
    # Aenderung der rechten Seite durch die Randbedingungen ------------- #    
    for i = 1:nc        
        # Boundaries         
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        # West boundary
        rhs[i]  +=  4*a*BC.val.W * DirW - 
                        2*a*BC.val.W*Δx * NeuW +
                        4*a*BC.val.E * DirE + 
                        2*a*BC.val.E*Δx * NeuE
    end
    # ------------------------------------------------------------------- #    
    # Compute new temperature ------------------------------------------- #
    cna.T      .=    K1 \ rhs
    # ------------------------------------------------------------------- #
end
# ======================================================================= #
# Time-dependent solvers, variable thermal parameters =================== #
# ======================================================================= #
"""
    ForwardEuler1D!(T, Py, Δt, Δy, nc, BC)

Solves the one-dimensional transient heat equation using an explicit Forward
Euler finite-difference scheme with variable thermal parameters and internal
heat production.

Temperature is defined at cell centroids, while thermal conductivity is
defined at the cell boundaries. The conductive heat flux is therefore
discretized in conservative form. Density, specific heat capacity, and heat
production are defined at the cell centroids.

Dirichlet and Neumann boundary conditions are supported at the southern and
northern boundaries. Ghost nodes are used to impose the boundary conditions.

# Arguments

    T       : Structure or tuple containing:
              `T`, the temperature array defined at the cell centroids, and
              `T_ex`, the extended temperature array including ghost nodes.
              Both arrays are updated in place.
    Py      : Structure or tuple containing the thermal parameters:
              `k`, thermal conductivity defined at the cell boundaries;
              `ρ`, density defined at the cell centroids;
              `cp`, specific heat capacity defined at the cell centroids;
              `H`, specific heat-production rate defined at the cell
              centroids.
    Δt      : Time step.
    Δy      : Grid spacing.
    nc      : Number of centroid nodes.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the southern and northern boundaries.

# Notes

Scalar values may be supplied for `k`, `ρ`, `cp`, and `H`. In this case, the
values are expanded internally to arrays of the required size.

The heat-production parameter `H` is expressed per unit mass. The corresponding
volumetric heat-production rate is therefore

    Q = ρ * H.

The method is first-order accurate in time and second-order accurate in space.
Because the temperature update is explicit, the timestep must satisfy the
diffusive stability condition. For variable thermal parameters, the timestep
should be selected using the largest local thermal diffusivity.

After each timestep, the updated centroid temperatures are copied to the
interior of `T.T_ex`. The ghost-node values are updated again when the routine
is called for the next timestep.

"""
function ForwardEuler1D!(T,Py,Δt,Δy,nc,BC)
    
    if size(Py.k,1) == 1
        k   =   Py.k.*ones(nc+1,1)
        ρ   =   Py.ρ.*ones(nc,1)
        cp  =   Py.cp.*ones(nc,1)
    else
        k   =   Py.k
        ρ   =   Py.ρ
        cp  =   Py.cp
    end
    if size(Py.H,1) == 1
        H   =   Py.H.*ones(nc,1)      #   [H] = W/kg; [Q] = [ρ*H], [Q] = W/m³
    else
        H   =   Py.H
    end

    # Define boundary conditions ---------------------------------------- #
    # South ---
    T.T_ex[1]   =   (BC.type.S==:Dirichlet) * (2 * BC.val.S - T.T_ex[2]) + 
                    (BC.type.S==:Neumann) * (T.T_ex[2] - BC.val.S*Δy)
    # North ---
    T.T_ex[end] =   (BC.type.N==:Dirichlet) * (2 * BC.val.N - T.T_ex[nc+1]) +
                    (BC.type.N==:Neumann) * (T.T_ex[nc+1] + BC.val.N*Δy)
    
    for j = 1:nc
        a       =   Δt/(Δy^2*ρ[j]*cp[j])
        T.T[j]  =   a*k[j]*T.T_ex[j] + 
                    (1-a*(k[j+1]+k[j]))*T.T_ex[j+1] +
                    a*k[j+1]*T.T_ex[j+2] +
                    H[j]*Δt/cp[j]
    end
    T.T_ex[2:end-1]     .= T.T
end

"""
    ComputeResiduals1D!(
        R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δx, Δt;
        C=0.0
    )

Computes the residual of the one-dimensional transient heat equation with
variable thermal properties and volumetric heat production.

The conductive heat flux is discretized in conservative form, with temperature
defined at the cell centroids and thermal conductivity at the cell boundaries.
Ghost nodes are used to impose Dirichlet and Neumann boundary conditions.

The routine evaluates the residual for the generalized θ-method, allowing
Backward Euler, Crank–Nicolson, and Forward Euler time discretizations through
the weighting parameter `C`. The residual is intended for defect-correction
iterations and is used together with the coefficient matrix assembled by
`AssembleMatrix1D!`.

# Arguments

    R       : Residual vector defined at the cell centroids.
    T       : Temperature at the current iteration or new time level.
    T_ex    : Extended current-temperature array including ghost nodes.
    T0      : Temperature at the previous time level.
    T_ex0   : Extended previous-temperature array including ghost nodes.
    Q       : Volumetric heat-production rate defined at the cell centroids.
    ∂T      : Structure or tuple containing the temperature gradients
              `∂x` and `∂x0`.
    q       : Structure or tuple containing the conductive heat fluxes
              `x` and `x0`.
    ρ       : Density defined at the cell centroids.
    Cp      : Specific heat capacity defined at the cell centroids.
    k       : Thermal conductivity defined at the cell boundaries.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western and eastern boundaries.
    Δx      : Grid spacing.
    Δt      : Time step.

# Keyword Arguments

    C       : Temporal weighting parameter (default: `0.0`):
              `0.0` for Backward Euler,
              `0.5` for Crank–Nicolson,
              `1.0` for Forward Euler.

# Notes

The residual is evaluated in conservative form as

- the transient storage term,
- the divergence of the conductive heat flux, and
- the volumetric heat-production term.

Backward Euler evaluates the conductive heat flux entirely at the new time
level, Crank–Nicolson averages the previous and current time levels, and
Forward Euler uses only the previous time level.

The routine updates the ghost-node values, temperature gradients, and heat
fluxes before evaluating the residual.
    
"""
function ComputeResiduals1D!(R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δx, Δt;C=0)
    if C < 1
        @. T_ex[2:end-1]    = T 
        T_ex[1]          = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2]) + (BC.type.W==:Neumann) * (T_ex[    2] - Δx/k[  1]*BC.val.W)
        T_ex[end]        = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1] + Δx/k[end]*BC.val.E)
        @. ∂T.∂x = (T_ex[2:end] - T_ex[1:end-1])/Δx
        @. q.x   = -k * ∂T.∂x
        if C==0.5
            @. T_ex0[2:end-1]   = T0 
            T_ex0[1]         = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2]) + (BC.type.W==:Neumann) * (T_ex0[    2] - Δx/k[  1]*BC.val.W)
            T_ex0[end]       = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1] + Δx/k[end]*BC.val.E)
            @. ∂T.∂x0   =   (T_ex0[2:end] - T_ex0[1:end-1])/Δx
            @. q.x0     =   -k * ∂T.∂x0
        end
    else
        @. T_ex0[2:end-1]   = T0 
        T_ex0[1]         = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2]) + (BC.type.W==:Neumann) * (T_ex0[    2] - Δx/k[  1]*BC.val.W)
        T_ex0[end]       = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1] + Δx/k[end]*BC.val.E)
        @. ∂T.∂x0   =   (T_ex0[2:end] - T_ex0[1:end-1])/Δx
        @. q.x0     =   -k * ∂T.∂x0
        @. q.x      =   q.x0
    end
    @. R     = ρ*Cp*(T - T0)/Δt + 
                    (1-C)*((q.x[2:end] - q.x[1:end-1])/Δx) + 
                    C*((q.x0[2:end] - q.x0[1:end-1])/Δx) - 
                    Q
end

"""
    AssembleMatrix1D(ρ, cp, k, Δx, Δt, nc, BC, K; C=0.0)

Assembles the coefficient matrix for the one-dimensional transient heat
equation with variable thermal properties.

Temperature is defined at the cell centroids, while thermal conductivity is
defined at the cell boundaries. The conductive term is discretized in
conservative flux-divergence form using second-order central finite
differences. Dirichlet and Neumann boundary conditions are incorporated
directly into the matrix coefficients at the western and eastern boundaries.

The temporal discretization is controlled by `C`. For `C < 1`, the matrix
contains the implicit contribution of the conductive term. For `C = 1`, the
conductive contribution vanishes and the matrix contains only the transient
storage term.

# Arguments

    ρ       : Density defined at the cell centroids.
    cp      : Specific heat capacity defined at the cell centroids.
    k       : Thermal conductivity defined at the cell boundaries.
    Δx      : Grid spacing.
    Δt      : Time step.
    nc      : Number of centroid nodes.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western and eastern boundaries.
    K       : Coefficient matrix to be assembled.

# Keyword Arguments

    C       : Temporal weighting parameter (default: `0.0`):
              `0.0` for Backward Euler,
              `0.5` for Crank–Nicolson,
              `1.0` for Forward Euler.

# Notes

The matrix corresponds to the current-time contribution of the generalized
time-discretization scheme used by `ComputeResiduals1D!`.

The assembled matrix is tridiagonal because each centroid temperature is
coupled only to its western and eastern neighbours. Boundary conditions modify
the coefficients in the first and last matrix rows.

The matrix is modified in place and finalized using `flush!(K)`.
"""
function AssembleMatrix1D( ρ, cp, k, Δx, Δt, nc, BC, K;C=0 )
    for i=1:nc
        # Equation number
        ii = i
        # Stencil
        iW = ii - 1
        iC = ii
        iE = ii + 1
        # Boundaries
        inW    =  i==1    ? false  : true   
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        inE    =  i==nc ? false  : true   
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        # Material coefficient
        kC = k[i]*(1-C)
        kE = k[i+1]*(1-C)
        # Linear system coefficients
        if inW K[ii,iW] = kC .* (DirW + NeuW - 1) ./ Δx .^ 2 end
        K[ii,iC] = cp[i] .* ρ[i] ./ Δt + (-kE .* (-DirE + NeuE - 1) ./ Δx + kC .* (DirW - NeuW + 1) ./ Δx) ./ Δx
        if inE K[ii,iE] = -kE .* (-DirE - NeuE + 1) ./ Δx .^ 2 end
    end
    return flush!(K)
end


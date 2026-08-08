using Base.Threads
using GeoModBox

"""
    IniTemperature!(type,M,NC,D,x,y;Tb = 600.0,Ta = 1200.0,σ  = 0.1)

Initialize a two-dimensional temperature field on the extended
cell-centered finite-difference grid.

The selected initial condition is assigned to `D.T_ex`, including its
ghost-cell layers. The interior temperature field `D.T` is subsequently
updated from

```julia
D.T_ex[2:end-1, 2:end-1]
```

If available, the auxiliary fields `D.T0`, `D.Told_ex`, `D.T_exD0`, and
`D.T_ex0` are initialized from the resulting temperature field.

# Arguments

    - `type`: Symbol selecting the initial temperature distribution.
    - `M`   : Model geometry containing `xmin`, `xmax`, `ymin`, and `ymax`.
    - `NC`  : Number of cell centers in the horizontal and vertical directions.
    - `D`   : Data structure containing the temperature fields.
    - `x`   : Horizontal coordinates, including the extended centroid coordinates
                `x.ce`.
    - `y`   : Vertical coordinates, including the extended centroid coordinates
                `y.ce`.

# Keyword arguments

    - `Tb=600.0`    : Background temperature for localized anomalies or bottom
                        temperature for vertically varying profiles.
    - `Ta=1200.0`   : Anomaly temperature for localized anomalies or top
                        temperature for vertically varying profiles.
    - `σ=0.1`       : Nondimensional width of the Gaussian anomaly.

# Initial-condition types

    - `:const`      : Spatially uniform temperature equal to `Tb`.
    - `:circle`     : Circular or elliptical region with temperature `Ta` embedded
                        in a background temperature `Tb`.
    - `:gaussian`   : Gaussian temperature anomaly with background temperature
                        `Tb` and peak temperature `Ta`.
    - `:block`      : Rectangular region with temperature `Ta` embedded in a
                        background temperature `Tb`.
    - `:linear`     : Linear conductive temperature profile between the prescribed
                        bottom and top temperatures.
    - `:lineara`    : Linear conductive temperature profile with an additional
                        elliptical temperature anomaly.
    - `:blankenbach`: Conductive temperature profile with the small sinusoidal
                        perturbation used to initialize the Blankenbach convection benchmark.

# Example

```julia
IniTemperature!(:circle,M,NC,D,x,y;Tb = 1200.0,Ta = 0.0,)
```
"""
@views function IniTemperature!(type,M,NC,D,x,y;Tb=600.0,Ta=1200.0,σ=0.1)
    if type==:const 
        # Constant temperature ---
        @. D.T_ex = Tb
    elseif type==:circle 
        # Circle shaped anomaly ---
        # Bereich der Anomalie ---       
        ri          =   .2
        xc          =   (M.xmax-M.xmin)/4
        yc          =   (M.ymax-M.ymin)/2
        α           =   0.0
        a_ell       =   .15*(M.ymax-M.ymin)
        b_ell       =   .15*(M.ymax-M.ymin)
        @threads for i = 1:NC.x+2 
            for j = 1:NC.y+2
                x_ell   =  x.ce[i]*cosd(α) + y.ce[j]*sind(α)
                y_ell   =  -x.ce[i]*sind(α) + y.ce[j]*cosd(α)
                Elli    =   ((x_ell - xc)/ a_ell)^2 + ((y_ell-yc)/ b_ell)^2
                if Elli <= ri 
                    D.T_ex[i,j]    =   Ta
                else
                    D.T_ex[i,j]    =   Tb
                end
            end
        end    
    elseif type==:gaussian        
        # κ           =   1e-6
        # AnalyticalSolution2D!(D.T, x.c, y.c, 0.0, (T0=Ampl,K=κ,σ=σ))
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                D.T_ex[i,j]    =   Tb + (Ta - Tb)*exp(-((x.ce[i]/((M.ymax-M.ymin)) - 0.20)^2 + 
                                    (y.ce[j]/((M.ymax-M.ymin)) - 0.5)^2)/σ^2)
            end
        end        
    elseif type==:block        
        # Bereich der Temperatur Anomalie ---
        xTl     =   M.xmin + (M.xmax-M.xmin)/8.0
        xTr     =   xTl + (M.xmax-M.xmin)/10.0
        yTu     =   M.ymin + (M.ymax-M.ymin)/2.0 - (M.ymax-M.ymin)/10.0 
        yTo     =   M.ymin + (M.ymax-M.ymin)/2.0 + (M.ymax-M.ymin)/10.0 
        # Anfangstemperatur Verteilung ---
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                if y.ce[j]>=yTu && y.ce[j] <= yTo && x.ce[i]>=xTl && x.ce[i]<=xTr
                    D.T_ex[i,j]    =   Ta
                else
                    D.T_ex[i,j]    =   Tb
                end
            end
        end        
    elseif type==:linear
        Ttop    =   Ta
        Tbot    =   Tb
        Tgrad   =   (Tbot-Ttop)/(M.ymax-M.ymin)         # [ K/m ]
        # @show Tgrad
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                D.T_ex[i,j] = -Tgrad*(y.ce[j]) + Ttop
            end
        end
    elseif type==:lineara
        # Bereich der Anomalie ---       
        ri          =   .1
        xc          =   (M.xmax-M.xmin)/4
        yc          =   -(M.ymax-M.ymin)/2
        α           =   0.0
        a_ell       =   3.0*ri*(M.ymax-M.ymin)
        b_ell       =   ri*(M.ymax-M.ymin)
        # --- 
        Ttop    =   Ta
        Tbot    =   Tb
        Tgrad   =   (Tbot-Ttop)/(M.ymax-M.ymin)         # [ K/m ]
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                x_ell   =  (x.ce[i]-xc)*cosd(α) + (y.ce[j]-yc)*sind(α)
                y_ell   =  -(x.ce[i]-xc)*sind(α) + (y.ce[j]-yc)*cosd(α)
                Elli    =   ((x_ell)/ a_ell)^2 + ((y_ell)/ b_ell)^2
                D.T_ex[i,j] = -Tgrad*(y.ce[j]) + Ttop
                if Elli <= 1.0
                    D.T_ex[i,j]    +=   0.2*D.T_ex[i,j]
                end
            end
        end
    elseif type == :blankenbach
        Ttop = Ta
        Tbot = Tb
        H    = M.ymax - M.ymin
        L    = M.xmax - M.xmin
        A    = 0.01
        # ---
        @threads for i = 1:NC.x+2
            for j = 1:NC.y+2
                xn = (x.ce[i] - M.xmin) / L
                zn = (y.ce[j] - M.ymin) / H

                Tlinear = Tbot + (Ttop - Tbot) * zn

                D.T_ex[i,j] =
                    Tlinear +
                    A * (Tbot - Ttop) *
                    cos(π*xn) * sin(π*zn)
            end
        end
    else
        throw(ArgumentError(
            "Unsupported initial temperature type `$type`. " *
            "Valid types are :const, :circle, :gaussian, :block, " *
            ":linear, :lineara, and :blankenbach."
        ))
    end
    # Assign temperature to regular field ---
    D.T         .=  D.T_ex[2:end-1,2:end-1]
    if  hasproperty(D, :T0)
        @. D.T0     =   D.T
    end
    if hasproperty(D, :Told_ex)
        @. D.Told_ex    =   D.T_ex
    end
    if hasproperty(D, :T_exD0)
        @. D.T_exD0 =   D.T_ex
    end
    if hasproperty(D, :T_ex0)
        @. D.T_ex0  =   D.T_ex
    end
    return D
end    

"""
    IniVelocity!(type,D,BC,NV,Δ,M,x,y;ε = 1e-15)

Initialize a predefined two-dimensional velocity field or update the
velocity boundary values on a staggered finite-difference grid.

For the rigid-body and shear-cell configurations, the complete staggered
velocity fields `D.vx` and `D.vy` are initialized. For simple-shear and
pure-shear configurations, the corresponding entries of `BC.val` are
updated and written to the velocity nodes and ghost nodes located along the
model boundaries. Interior velocity values are retained in these cases,
allowing the function to update deformation boundary conditions after a
grid-remeshing step without erasing the current Stokes solution.

# Arguments

    - `type`: Symbol selecting the velocity configuration.
    - `D`   : Data structure containing the staggered velocity fields `vx` and
                `vy`.
    - `BC`  : Velocity boundary-condition structure containing `type` and `val`.
    - `NV`  : Number of grid vertices in the horizontal and vertical directions.
    - `Δ`   : Grid spacing containing `x` and `y`.
    - `M`   : Model geometry containing `xmin`, `xmax`, `ymin`, and `ymax`.
    - `x`   : Horizontal coordinates at cell centers, vertices, and staggered
                velocity locations.
    - `y`   : Vertical coordinates at cell centers, vertices, and staggered
                velocity locations.

# Keyword arguments

    - `ε=1e-15` : Background deformation rate in s⁻¹ used by the simple- and
                    pure-shear velocity configurations.

# Velocity configurations

    - `:RigidBody`  : Clockwise rigid-body rotation within a circular region.
                        The velocity decreases to zero outside the prescribed rotating region.
    - `:ShearCell`  : Analytic cellular velocity field used for advection tests.
    - `:SimpleShear`: Simple-shear deformation with horizontal velocity varying
                        linearly with vertical position and zero vertical velocity.
    - `:PureShear`  : Domain-centered incompressible pure shear,

  ```math
  v_x = -\\dot{\\varepsilon}(x-x_c), \\qquad
  v_y =  \\dot{\\varepsilon}(y-y_c).
  ```

    - `:ShearBandPS`: Origin-centered pure-shear field used by the deforming
                        shear-band experiment,

  ```math
  v_x = -\\dot{\\varepsilon}x, \\qquad
  v_y =  \\dot{\\varepsilon}y.
  ```

For `:SimpleShear`, `:PureShear`, and `:ShearBandPS`, the boundary types in
`BC.type` must be configured consistently with the selected deformation
field before calling this function.

# Returns

Returns the modified data and boundary-condition structures as `(D, BC)`.

# Example

```julia
D, VBC = IniVelocity!(:PureShear,D,VBC,NV,Δ,M,x,y;ε = 1e-15)
```
"""
@views function IniVelocity!(type,D,BC,NV,Δ,M,x,y;ε=1e-15)
    xc  =   0.5 * (M.xmin + M.xmax)
    yc  =   0.5 * (M.ymin + M.ymax)
    L   =   M.xmax - M.xmin
    H   =   M.ymax - M.ymin

    if type==:RigidBody
        # Rigid Body Rotation ---
        # We assume a maximum and minimum velocity of 0.5 cm/a, respectively! 
        @threads for i = 1:NV.x
            for j = 1:NV.y+1
                D.vx[i,j]  =    (y.ce[j]-yc) / H
            end
        end
        @threads for i = 1:NV.x+1
            for j = 1:NV.y
                D.vy[i,j]  =   -(x.ce[i]-xc) / H
            end
        end
        
        Radx        =   zeros(size(D.vx))
        Rady        =   zeros(size(D.vy))

        @. Radx     =   sqrt((x.vx2d - xc)^2 + (y.vx2d - yc)^2)
        @. Rady     =   sqrt((x.vy2d - xc)^2 + (y.vy2d - yc)^2)

        R           =   0.5 * min(L, H) - max(Δ.x, Δ.y)

        D.vx[Radx .> R] .= 0.0
        D.vy[Rady .> R] .= 0.0

        @. D.vx     =   D.vx/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
        @. D.vy     =   D.vy/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
        
    elseif type==:ShearCell
        # Convection Cell with a Shear Deformation --- (REF?!)
        @threads for i in 1:NV.x
            for j in 1:NV.y+1
                xn = (x.v[i]  - M.xmin) / L
                yn = (y.ce[j] - M.ymin) / H

                D.vx[i,j] = -(L/H) * sinpi(xn) * cospi(yn)
            end
        end

        @threads for i in 1:NV.x+1
            for j in 1:NV.y
                xn = (x.ce[i] - M.xmin) / L
                yn = (y.v[j]  - M.ymin) / H

                D.vy[i,j] = cospi(xn) * sinpi(yn)
            end
        end
        @. D.vx     =   D.vx/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
        @. D.vy     =   D.vy/(100.0*(60.0*60.0*24.0*365.25))        # [m/s]
    elseif type==:SimpleShear || type==:PureShear || type==:ShearBandPS
        if type==:SimpleShear
            # Horizontal velocity 
            @. BC.val.S     =   (y.v2d[:,1]   - M.ymin) * ε     #   South
            @. BC.val.N     =   (y.v2d[:,end] - M.ymin) * ε     #   North
            @. BC.val.vxW   =   (y.c2d[1,:]   - M.ymin) * ε     #   West    
            @. BC.val.vxE   =   (y.c2d[end,:] - M.ymin) * ε     #   East
            # Vertical velocity 
            @. BC.val.vyS   =   0.0
            @. BC.val.vyN   =   0.0
            @. BC.val.W     =   0.0
            @. BC.val.E     =   0.0            
            
        elseif type==:PureShear
            # Tangential vx values at South and North
            @. BC.val.S     =   -(x.vx2d[:,1]   - xc) * ε
            @. BC.val.N     =   -(x.vx2d[:,end] - xc) * ε

            # Normal vx values at West and East
            @. BC.val.vxW   =   -(x.vx2d[1,2:end-1]   - xc) * ε
            @. BC.val.vxE   =   -(x.vx2d[end,2:end-1] - xc) * ε

            # Normal vy values at South and North
            @. BC.val.vyS   =   (y.vy2d[2:end-1,1]   - yc) * ε
            @. BC.val.vyN   =   (y.vy2d[2:end-1,end] - yc) * ε

            # Tangential vy values at West and East
            @. BC.val.W     =   (y.v2d[1,:]   - yc) * ε
            @. BC.val.E     =   (y.v2d[end,:] - yc) * ε

            # # Horizontal velocity 
            # @. BC.val.S     =   -(x.vx2d[:,1]-(M.xmax-M.xmin)/2)*ε               #   South
            # @. BC.val.N     =   -(x.vx2d[:,end]-(M.xmax-M.xmin)/2)*ε             #   North
            # @. BC.val.vxW   =   -(x.vx2d[1,2:end-1]-(M.xmax-M.xmin)/2)*ε         #   West
            # @. BC.val.vxE   =   -(x.vx2d[end,2:end-1]-(M.xmax-M.xmin)/2)*ε       #   East
            # # Vertical velocity 
            # @. BC.val.vyS   =   (y.vy2d[2:end-1,1]+(M.ymax-M.ymin)/2)*ε         #   South
            # @. BC.val.vyN   =   (y.vy2d[2:end-1,end]+(M.ymax-M.ymin)/2)*ε       #   North
            # @. BC.val.W     =   (y.v2d[1,:]+(M.ymax-M.ymin)/2)*ε                #   West
            # @. BC.val.E     =   (y.v2d[end,:]+(M.ymax-M.ymin)/2)*ε              #   East
        elseif type==:ShearBandPS
            # Horizontal velocity 
            @. BC.val.S     =   -x.vx2d[:,1]*ε                                  #   South
            @. BC.val.N     =   -x.vx2d[:,end]*ε                                #   North
            @. BC.val.vxW   =   -x.vx2d[1,2:end-1]*ε                            #   West
            @. BC.val.vxE   =   -x.vx2d[end,2:end-1]*ε                          #   East
            # Vertical velocity 
            @. BC.val.vyS   =   y.vy2d[2:end-1,1]*ε                             #   South
            @. BC.val.vyN   =   y.vy2d[2:end-1,end]*ε                           #   North
            @. BC.val.W     =   y.v2d[1,:]*ε                                    #   West
            @. BC.val.E     =   y.v2d[end,:]*ε                                  #   East
        end
        D.vx[1,2:end-1]     .=   BC.val.vxW     #   West - vx
        D.vy[1,:]           .=   BC.val.W       #   West - vy
        D.vx[end,2:end-1]   .=   BC.val.vxE     #   East - vx
        D.vy[end,:]         .=   BC.val.E       #   East - vy
        D.vy[2:end-1,1]     .=   BC.val.vyS     #   South - vy
        D.vx[:,1]           .=   BC.val.S       #   South - vx
        D.vy[2:end-1,end]   .=   BC.val.vyN     #   North - vy
        D.vx[:,end]         .=   BC.val.N       #   North - vx
    else
        throw(ArgumentError(
            "Unsupported initial velocity type `$type`. " *
            "Valid types are :RigidBody, :ShearCell, :SimpleShear, " *
            ":PureShear, and :ShearBandPS."
        ))
    end
    return D, BC
end

"""
    IniPhase!(type, D, M, x, y, NC; phase=(0, 1))

Initialize a two-dimensional phase field at the cell centers of a regular
finite-difference grid.

!!! warning
    This function defines phases directly on the Eulerian grid and is
    retained primarily for compatibility with older examples. For models
    involving material transport, phases should preferably be initialized
    on tracers and subsequently interpolated to the grid.

# Arguments

    - `type`: Symbol selecting the initial phase distribution.
    - `D`   : Data structure containing the cell-centered phase field `p` and the
                extended phase field `p_ex`.
    - `M`   : Model geometry containing `xmin`, `xmax`, `ymin`, and `ymax`.
    - `x`   : Horizontal cell-center coordinates.
    - `y`   : Vertical cell-center coordinates.
    - `NC`  : Number of cell centers in the horizontal and vertical directions.

# Keyword arguments

    - `phase=(0, 1)`: Two phase values given as
                        `(background_phase, anomaly_phase)`.

# Supported phase distributions

    - `:block`: Rectangular anomaly located between 40% and 60% of the model
                    width and between 10% and 30% of the model height below the upper
                    boundary.

The interior phase field is assigned to `D.p`. The corresponding extended
field `D.p_ex` is then updated using constant extrapolation into the ghost
cells.

# Example

```julia
IniPhase!(:block,D,M,x,y,NC;phase = (0, 1))
"""
@views function IniPhase!(type,D,M,x,y,NC;phase=(0, 1))

    if type==:block
        length(phase) == 2 ||
        throw(ArgumentError(
            "`phase` must contain exactly two values: " *
            "`(background_phase, anomaly_phase)`."
        ))
        # Bereich der Anomalie ---
        L       =   M.xmax - M.xmin
        H       =   M.ymax - M.ymin
        xL      =   M.xmin + 2/5 * L
        xR      =   M.xmin + 3/5 * L
        yO      =   M.ymax - 0.1 * H
        yU      =   M.ymax - 0.3 * H        
        
        # Phase ---
        for i = 1:NC.x
            for j = 1:NC.y
                if y.c[j]>=yU && y.c[j] <= yO && x.c[i]>=xL && x.c[i]<=xR
                    D.p[i,j]    =   phase[2]    #   anomaly 
                else
                    D.p[i,j]    =   phase[1]    #   background
                end
            end
        end
        D.p_ex[2:end-1,2:end-1]     .=  D.p
        D.p_ex[1,:]     .=   D.p_ex[2,:]
        D.p_ex[end,:]   .=   D.p_ex[end-1,:]
        D.p_ex[:,1]     .=   D.p_ex[:,2]
        D.p_ex[:,end]   .=   D.p_ex[:,end-1]
    else
        throw(ArgumentError(
            "Unsupported initial phase type `$type`. " *
            "Currently supported type: :block."
        ))
    end
    return D
end
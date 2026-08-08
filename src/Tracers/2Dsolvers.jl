using Base.Threads
using Statistics, Printf

"""
    Ma = TMarkers(x, y, T, phase)

Container storing the Lagrangian marker properties used by the advection
routines.

Each marker stores its Cartesian coordinates, temperature, and phase ID.
The structure is primarily used for marker-based temperature advection and
for tracking material interfaces during thermo-mechanical simulations.

# Fields

    - `x`     : Horizontal marker coordinates.
    - `y`     : Vertical marker coordinates.
    - `T`     : Marker temperature.
    - `phase` : Integer phase identifier carried by each marker.

All arrays must have identical length corresponding to the total number of
markers.

"""
mutable struct TMarkers
    x       ::  Array{Float64,1}
    y       ::  Array{Float64,1}
    T       ::  Array{Float64,1}
    phase   ::  Array{Int64,1}
end

"""
    Ma = Markers(x, y, phase)

Container storing the Lagrangian marker positions and phase IDs.

The structure is primarily used for marker-based material advection. Each
marker stores its Cartesian coordinates together with an integer phase ID,
allowing the phase distribution to evolve independently of the Eulerian
computational grid.

# Fields

    - `x`     : Horizontal marker coordinates.
    - `y`     : Vertical marker coordinates.
    - `phase` : Integer phase identifier carried by each marker.

All arrays must have identical length corresponding to the total number of
markers.
"""
mutable struct Markers
    x       ::  Array{Float64,1}
    y       ::  Array{Float64,1}
    phase   ::  Array{Int64,1}
end

"""
    IniTracer2D(
        Aparam,
        nmx,
        nmy,
        Δ,
        M,
        NC,
        noise,
        ini,
        phase;
        xc,
        yc,
        λ,
        δA,
        ellA,
        ellB,
        α,
    )

Initialize a regular two-dimensional distribution of Lagrangian tracers and
assign their initial phase IDs.

The tracers are distributed uniformly within each finite-difference cell.
Depending on `Aparam`, the returned marker structure stores either phase
information only or both phase and temperature. An optional random
perturbation can be applied to the initial tracer coordinates.

# Arguments

    - `Aparam`  : Defines the marker properties:
                    - `:phase` initializes markers carrying phase IDs.
                    - `:thermal` initializes markers carrying phase IDs and temperature.
    - `nmx`     : Number of tracers per cell in the horizontal direction.
    - `nmy`     : Number of tracers per cell in the vertical direction.
    - `Δ`       : Grid spacing containing the horizontal and vertical increments.
    - `M`       : Model geometry containing the domain limits.
    - `NC`      : Number of finite-difference cells in the horizontal and vertical
                    directions.
    - `noise`   : Controls whether random perturbations are added to the tracer
                    positions. Use `1` to add noise and `0` to retain the regular distribution.
    - `ini`     : Symbol selecting the initial phase distribution.
    - `phase`   : Collection containing the phase IDs. `phase[1]` represents the
                    background or matrix phase, while `phase[2]` represents the anomaly,
                    inclusion, or second layer.

# Keyword arguments

    - `xc=(M.xmax-M.xmin)`  : Horizontal reference coordinate used by the
                                Rayleigh–Taylor instability setup.
    - `yc=(M.ymax-M.ymin)/2`: Vertical reference coordinate used by the
                                Rayleigh–Taylor instability setup.
    - `λ=1.0e3`             : Wavelength of the cosine perturbation used for `:RTI`.
    - `δA=5e2/15`           : Amplitude of the perturbation used for `:RTI`.
    - `ellA=100.0`          : Major semiaxis of the elliptical inclusion or radius of the
                                circular shear-band inclusion.
    - `ellB=100.0`          : Minor semiaxis of the elliptical inclusion.
    - `α=0.0`               : Rotation angle of the elliptical inclusion in degrees.

# Initial phase distributions

    - `:block`              : Rectangular anomaly embedded in a homogeneous background phase.
    - `:RTI`                : Two-layer configuration with a cosine perturbation of the tracer
                                positions for a Rayleigh–Taylor instability.
    - `:Inclusion`          : Elliptical inclusion located at the center of the model
                                domain. The inclusion may be rotated by the angle `α`.
    - `:ShearBandSetting`   : Circular weak inclusion centered at the lower
                                boundary for the thermo-mechanical shear-localization experiment.

# Returns

Returns either a [`Markers`](@ref) or [`TMarkers`](@ref) structure containing
the initialized tracer coordinates, phase IDs, and, where applicable,
temperature values.

# Example

```julia
Ma = IniTracer2D(
    :phase,
    nmx,
    nmy,
    Δ,
    M,
    NC,
    1,
    :RTI,
    [0, 1],
)
```
"""
@views function IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,ini,phase;
                    xc=(M.xmax-M.xmin),yc=(M.ymax-M.ymin)/2,λ=1.0e3,δA=5e2/15,
                    ellA=100.0,ellB=100.0,α=0.0)
    
    nmark   =   nmx*nmy*NC.x*NC.y

    # Initialise markers ---
    Δxm, Δym  =   Δ.x/nmx, Δ.y/nmy 
    
    xm  =   LinRange(M.xmin+Δxm/2.0, M.xmax-Δxm/2.0, NC.x*nmx)
    ym  =   LinRange(M.ymin+Δym/2.0, M.ymax-Δym/2.0, NC.y*nmy)        
    
    (xmi,ymi) = ([x for x=xm,y=ym], [y for x=xm,y=ym])
    
    # Over allocate markers ---
    if Aparam==:thermal        
        Ma  =   TMarkers( vec(xmi), vec(ymi), zeros(Float64, nmark), zeros(Int64, nmark) )
    elseif Aparam==:phase
        Ma  =   Markers( vec(xmi), vec(ymi), zeros(Int64, nmark) )
    end

    # add noise ---
    if noise==1
        @threads for k=1:nmark
            Ma.x[k] += (rand()-0.5)*Δxm # * 3.0 / 4.0
            Ma.y[k] += (rand()-0.5)*Δym # * 3.0 / 4.0
        end
    end

    if ini==:block
        # Geometry of the block anomaly ---
        xL      =   2/5 * (M.xmax-M.xmin)
        xR      =   3/5 * (M.xmax-M.xmin)
        yO      =   0.1 * (M.ymin-M.ymax)
        yU      =   0.3 * (M.ymin-M.ymax)        
        
        # phase ---
        for k = 1:nmark
            if Ma.y[k]>=yU && Ma.y[k]<=yO && Ma.x[k]>=xL && Ma.x[k]<=xR
                Ma.phase[k]     =   phase[2]    #   anomaly 
            else
                Ma.phase[k]     =   phase[1]    #   background
            end
        end
    elseif ini==:RTI
        @threads for k=1:nmark
            # Layer interface  --- 
            δAm     =   cos(2*π*((Ma.x[k] - 0.5*xc)/λ))*δA
            if abs(Ma.y[k]) >=  yc
                Ma.phase[k]     =   phase[2]    #   Lower layer
            else
                Ma.phase[k]     =   phase[1]    #   Upper layer
            end            
            Ma.y[k]     -= δAm
            if Ma.y[k] > M.ymax 
                Ma.y[k] += M.ymin
                Ma.phase[k] = phase[2]
            elseif Ma.y[k] < M.ymin
                Ma.y[k] -= M.ymin
                Ma.phase[k] = phase[1]
            end
        end
    elseif ini==:Inclusion
        # Circle shaped anomaly ---
        # Bereich der Anomalie ---       
        xc          =   (M.xmin+M.xmax)/2
        yc          =   (M.ymin+M.ymax)/2
        @threads for k = 1:nmark
            x_ell   =  (Ma.x[k]-xc)*cosd(α) + (Ma.y[k]-yc)*sind(α)
            y_ell   =  -(Ma.x[k]-xc)*sind(α) + (Ma.y[k]-yc)*cosd(α)
            Elli    =   (x_ell/ellA)^2 + (y_ell/ellB)^2
            if Elli <= 1
                Ma.phase[k]     =   phase[2]    #   Inclusion
            else
                Ma.phase[k]     =   phase[1]    #   Matrix
            end
        end
    elseif ini==:ShearBandSetting
        # Bereich der Anomalie ---       
        xc          =   (M.xmin+M.xmax)/2
        yc          =   M.ymin
        @threads for k = 1:nmark
            x_ell   =  Ma.x[k] - xc
            y_ell   =  Ma.y[k] - yc
            Elli    =   (x_ell/ellA)^2 + (y_ell/ellA)^2
            if Elli <= 1
                Ma.phase[k]     =   phase[2]    #   Inclusion
            else
                Ma.phase[k]     =   phase[1]    #   Matrix
            end
        end
    end

    return Ma
end

"""
    CountMPC(Ma, nmark, MPC, M, x, y, Δ, NC, NV; verbose=false)

Count the number of active markers associated with each cell center and
grid vertex.

Markers located outside the current model domain are deactivated by setting
their phase ID to `-1`. The remaining active markers are counted separately
for the cell-centered and vertex-centered grids. Thread-local counting
arrays are used during the parallel loops and are subsequently summed to
obtain the total marker counts stored in `MPC.c` and `MPC.v`.

# Arguments

- `Ma`      : Marker structure containing the marker coordinates and phase IDs.
- `nmark`   : Total number of markers.
- `MPC`     : Structure containing the marker-count arrays:
    - `MPC.c`   : Total number of active markers associated with each cell center.
    - `MPC.v`   : Total number of active markers associated with each grid vertex.
    - `MPC.th`  : Thread-local marker counts for the cell-centered grid.
    - `MPC.thv` : Thread-local marker counts for the vertex-centered grid.
- `M`       : Model geometry containing the domain limits.
- `x`       : Horizontal coordinates of the cell centers and vertices.
- `y`       : Vertical coordinates of the cell centers and vertices.
- `Δ`       : Grid spacing in the horizontal and vertical directions.
- `NC`      : Number of cell centers in each coordinate direction.
- `NV`      : Number of grid vertices in each coordinate direction.

# Keyword arguments

- `verbose=false`: Print the number of markers located outside the model
  domain.

# Notes

Only markers with

```julia
Ma.phase[k] >= 0
```
are included in the marker counts. Markers outside the domain are assigned
a phase ID of -1 and are excluded from subsequent interpolation and
advection operations.

The function modifies Ma.phase, MPC.c, MPC.v, MPC.th, and
MPC.thv in place.
```
# Example

CountMPC(
    Ma,
    nmark,
    MPC,
    M,
    x,
    y,
    Δ,
    NC,
    NV;
    verbose = false,
)
"""
@views function CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV;verbose=false)
     # Disable markers outside of the domain
     @threads for k=1:nmark
        if (Ma.x[k]<M.xmin || Ma.x[k]>M.xmax || Ma.y[k]<M.ymin || Ma.y[k]>M.ymax) 
            @inbounds Ma.phase[k] = -1
        end
    end
    # How many are outside? save indices for reuse  
    nmark_out_th    =   zeros(Int64, nthreads())  
    @threads for k=1:nmark
        if Ma.phase[k] == -1
            nmark_out_th[threadid()] += 1
        end
    end
    nmark_out       =   0
    for ith=1:nthreads()
        nmark_out += nmark_out_th[ith]
    end
    if verbose 
        @printf("%d markers out\n", nmark_out)
    end

    # Initialize marker per cell per thread array ---
    @threads for j = 1:NC.y
        for i = 1:NC.x
            for ith=1:nthreads()
                MPC.th[ith,i,j] = 0.0
            end
        end
    end

    # Count marker per cell per thread ---
    @threads for k=1:nmark
        if (Ma.phase[k]>=0)
            # Get the column:
            dstx = Ma.x[k] - x.c[1]
            i    = Int64(round(ceil( (dstx/Δ.x) + 0.5)))
            # Get the line:
            dsty = Ma.y[k] - y.c[1]
            j    = Int64(round(ceil( (dsty/Δ.y) + 0.5)))
            # Increment cell count
            MPC.th[threadid(),i,j] += 1.0
        end
    end

    @threads for j=1:NC.y
        for i=1:NC.x
            for ith=1:nthreads()
                if ith == 1 
                    MPC.c[i,j] = 0.0
                end
                MPC.c[i,j] += MPC.th[ith,i,j]
            end
        end
    end

    # Initialize marker per vertex per thread array ---
    @threads for j = 1:NV.y
        for i = 1:NV.x
            for ith=1:nthreads()
                MPC.thv[ith,i,j] = 0.0
            end
        end
    end

    # Count marker per vertex per thread ---
    @threads for k=1:nmark
        if (Ma.phase[k]>=0)
            # Get the column:
            dstx = Ma.x[k] - x.v[1]
            i    = Int64(round(ceil( (dstx/Δ.x) + 0.5)))
            # Get the line:
            dsty = Ma.y[k] - y.v[1]
            j    = Int64(round(ceil( (dsty/Δ.y) + 0.5)))
            # Increment cell count
            MPC.thv[threadid(),i,j] += 1.0
        end
    end

    @threads for j=1:NV.y
        for i=1:NV.x
            for ith=1:nthreads()
                if ith == 1 
                    MPC.v[i,j] = 0.0
                end
                MPC.v[i,j] += MPC.thv[ith,i,j]
            end
        end
    end

    # MPC.min[it]         = minimum(MPC.cv)
    # MPC.max[it]         = maximum(MPC.cv)
    # MPC.mean[it]        = mean(MPC.cv)
    #MPC.tot_reseed[it] = nmark_add
    # return Ma
end

"""
    VxFromVxNodes(Vx, k, Ma, x, y, Δ, NC, new)

Interpolate the horizontal velocity component from the staggered
`vₓ` nodes to the position of marker `k`.

The function first determines the surrounding staggered velocity nodes and
performs bilinear interpolation in the horizontal and vertical directions.
Marker coordinates located close to the model boundaries are restricted to
valid interpolation indices.

An optional higher-order velocity correction can be activated through
`new`. The implementation is adapted from the marker-advection routines in
`M2Dpt_Julia`.

# Arguments

- `Vx`: Horizontal velocity field defined at the staggered `vₓ` nodes.
- `k`: Index of the marker for which the velocity is evaluated.
- `Ma`: Marker structure containing the marker coordinates.
- `x`: Horizontal coordinates of the staggered velocity nodes.
- `y`: Vertical coordinates of the extended cell centers.
- `Δ`: Grid spacing in the horizontal and vertical directions.
- `NC`: Number of cell centers in each coordinate direction.
- `new`: Switch controlling the interpolation formulation:
  - `0`: Use bilinear interpolation.
  - `1`: Apply the additional higher-order velocity correction where
    sufficient neighboring nodes are available.

# Returns

Returns the interpolated horizontal marker velocity `vxm`.

# Example

```julia
vxm = VxFromVxNodes(
    D.vx,
    k,
    Ma,
    x,
    y,
    Δ,
    NC,
    0,
)
```

# Reference

The interpolation formulation is modified from the marker-advection
implementation in M2Dpt_Julia:

https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
"""
@views function VxFromVxNodes(Vx, k, Ma, x, y, Δ, NC, new)
    # Interpolate vx
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (Ma.x[k] -  x.v[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1))
    if i<1
        i = 1
    elseif i>NC.x
        i = NC.x
    end
    if j<1
        j = 1
    elseif j> NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   Ma.x[k] -  x.v[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute vx velocity for the top and bottom of the cell --------- #
    vxm13   =   Vx[i,j  ] * (1-Δxmj/Δ.x) + Vx[i+1,j  ]*Δxmj/Δ.x
    vxm24   =   Vx[i,j+1] * (1-Δxmj/Δ.x) + Vx[i+1,j+1]*Δxmj/Δ.x
    # ---------------------------------------------------------------- #
    if new==1 
        if dxmj/dx>=0.5
            if i<ncx
                vxm13 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i,j  ] - 2.0*Vx[i+1,j  ] + Vx[i+2,j  ])
                vxm24 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i,j+1] - 2.0*Vx[i+1,j+1] + Vx[i+2,j+1])
            end
        else
            if i>1
                vxm13 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i-1,j  ] - 2.0*Vx[i,j  ] + Vx[i+1,j  ])
                vxm24 += 0.5*((dxmj/Δ.x-0.5)^2) * (Vx[i-1,j+1] - 2.0*Vx[i,j+1] + Vx[i+1,j+1])
            end
        end
    end
    # ----------------------------------------------------------------- #
    # Compute vx ------------------------------------------------------ #
    vxm     =   (1-Δymi/Δ.y) * vxm13 + (Δymi/Δ.y) * vxm24
    # ----------------------------------------------------------------- #
    return vxm
end

"""
    VyFromVyNodes(Vy, k, Ma, x, y, Δ, NC, new)

Interpolate the vertical velocity component from the staggered
`vᵧ` nodes to the position of marker `k`.

The function first identifies the surrounding staggered velocity nodes and
computes the marker velocity using bilinear interpolation. Marker
coordinates close to the model boundaries are clamped to valid
interpolation indices to ensure that only existing grid nodes are used.

An optional higher-order interpolation correction can be enabled through
`new`. The implementation is adapted from the marker-advection routines in
`M2Dpt_Julia`.

# Arguments

- `Vy`: Vertical velocity field defined at the staggered `vᵧ` nodes.
- `k`: Index of the marker for which the velocity is evaluated.
- `Ma`: Marker structure containing the marker coordinates.
- `x`: Horizontal coordinates of the extended cell centers.
- `y`: Vertical coordinates of the staggered velocity nodes.
- `Δ`: Grid spacing in the horizontal and vertical directions.
- `NC`: Number of cell centers in each coordinate direction.
- `new`: Switch controlling the interpolation formulation:
  - `0`: Use bilinear interpolation.
  - `1`: Apply the additional higher-order velocity correction where
    sufficient neighboring nodes are available.

# Returns

Returns the interpolated vertical marker velocity `vym`.

# Example

```julia
vym = VyFromVyNodes(
    D.vy,
    k,
    Ma,
    x,
    y,
    Δ,
    NC,
    0,
)
```
# Reference

The interpolation formulation is modified from the marker-advection
implementation in M2Dpt_Julia:

https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
"""
@views function VyFromVyNodes(Vy, k, Ma, x, y, Δ, NC, new)
    # Interpolate vy
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (Ma.y[k] -  y.v[1])/Δ.y ) + 1))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j = 1
    elseif j>NC.y
        j = NC.y
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    dxmj = Ma.x[k] - x.ce[i]
    dymi = Ma.y[k] -  y.v[j]
    # ---------------------------------------------------------------- #
    # Compute vy velocity for the left and right of the cell --------- #
    vym12 = Vy[i,j  ]*(1-dymi/Δ.y) + Vy[i  ,j+1]*dymi/Δ.y
    vym34 = Vy[i+1,j]*(1-dymi/Δ.y) + Vy[i+1,j+1]*dymi/Δ.y
    # ---------------------------------------------------------------- #
    if new==1 
        if dymi/dy>=0.5
            if j<ncy
                vym12 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i,j  ] - 2.0*Vy[i,j+1  ] + Vy[i,j+2  ]);
                vym34 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i+1,j] - 2.0*Vy[i+1,j+1] + Vy[i+1,j+2]);
            end      
        else
            if j>1
                vym12 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i,j-1  ] - 2.0*Vy[i,j  ] + Vy[i,j+1  ]);
                vym34 += 0.5*((dymi/Δ.y-0.5)^2) * ( Vy[i+1,j-1] - 2.0*Vy[i+1,j] + Vy[i+1,j+1]);
            end
        end
    end
    # ---------------------------------------------------------------- #
    # Compute vy ----------------------------------------------------- #
    vym = (1-dxmj/Δ.x)*vym12 + (dxmj/Δ.x)*vym34
    # ---------------------------------------------------------------- #
    return vym
end

"""
    VxVyFromPrNodes(Vxp, Vyp, k, Ma, x, y, Δ, NC)

Interpolate the horizontal and vertical velocity components from the
cell-centered grid to the position of marker `k`.

The function identifies the four surrounding cell-center nodes and applies
bilinear interpolation to both velocity components. Marker coordinates
close to the model boundaries are restricted to valid interpolation
indices.

# Arguments

- `Vxp`: Horizontal velocity field defined at the cell centers.
- `Vyp`: Vertical velocity field defined at the cell centers.
- `k`: Index of the marker for which the velocity is evaluated.
- `Ma`: Marker structure containing the marker coordinates.
- `x`: Horizontal coordinates of the extended cell centers.
- `y`: Vertical coordinates of the extended cell centers.
- `Δ`: Grid spacing in the horizontal and vertical directions.
- `NC`: Number of cell centers in each coordinate direction.

# Returns

Returns the interpolated horizontal and vertical marker velocities as

```julia
vxm, vym
```

# Example 
vxm, vym = VxVyFromPrNodes(
    D.vxc,
    D.vyc,
    k,
    Ma,
    x,
    y,
    Δ,
    NC,
)

# Reference

The interpolation formulation is modified from the marker-advection
implementation in M2Dpt_Julia:

https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
"""
@views function VxVyFromPrNodes(Vxp ,Vyp, k, Ma, x, y, Δ, NC)
    # Interpolate vx, vy
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1.0))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j=1
    elseif j>NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   Ma.x[k] - x.ce[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute weights ------------------------------------------------ #
    wtmij   =   (1.0-Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j  =   (1.0-Δxmj/Δ.x)*(    Δymi/Δ.y)    
    wtmij1  =   (    Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j1 =   (    Δxmj/Δ.x)*(    Δymi/Δ.y)
    # ---------------------------------------------------------------- #
    # Compute vx, vy velocity ---------------------------------------- #
    vxm = Vxp[i,j]*wtmij + Vxp[i,j+1]*wtmi1j + Vxp[i+1,j]*wtmij1 + Vxp[i+1,j+1]*wtmi1j1
    vym = Vyp[i,j]*wtmij + Vyp[i,j+1]*wtmi1j + Vyp[i+1,j]*wtmij1 + Vyp[i+1,j+1]*wtmi1j1
    # ---------------------------------------------------------------- #
    return vxm, vym
end

"""
    FromCtoM(Prop, k, Ma, x, y, Δ, NC)

Interpolate a cell-centered property to the position of marker `k`.

The function identifies the four surrounding cell centers enclosing the
marker position and evaluates the property using bilinear interpolation.
Marker coordinates located close to the model boundaries are restricted to
valid interpolation indices to ensure that only existing grid values are
used.

This routine can be used to initialize or update marker properties from
Eulerian fields, such as temperature, density, viscosity, or any other
scalar quantity defined at the cell centers.

# Arguments

- `Prop`: Property field defined at the cell centers.
- `k`: Index of the marker for which the property is evaluated.
- `Ma`: Marker structure containing the marker coordinates.
- `x`: Horizontal coordinates of the extended cell centers.
- `y`: Vertical coordinates of the extended cell centers.
- `Δ`: Grid spacing in the horizontal and vertical directions.
- `NC`: Number of cell centers in each coordinate direction.

# Returns

Returns the interpolated property value at the position of marker `k`.

# Example

```julia
for k = 1:nmark
    Ma.T[k] = FromCtoM(D.T_ex, k, Ma, x, y, Δ, NC)
end
```
"""
@views function FromCtoM(Prop, k, Ma, x, y, Δ, NC)
    # Interpolate Property from Centroids to Marker ---
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1.0))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j=1
    elseif j>NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   Ma.x[k] - x.ce[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute weights ------------------------------------------------ #
    wtmij   =   (1.0-Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j  =   (1.0-Δxmj/Δ.x)*(    Δymi/Δ.y)    
    wtmij1  =   (    Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j1 =   (    Δxmj/Δ.x)*(    Δymi/Δ.y)
    # ---------------------------------------------------------------- #
    # Compute Marker Property ---------------------------------------- #
    Propm   =   Prop[i,j]*wtmij + 
                    Prop[i,j+1]*wtmi1j + 
                    Prop[i+1,j]*wtmij1 + 
                    Prop[i+1,j+1]*wtmi1j1
    # ---------------------------------------------------------------- #
    return Propm
end

"""
    Markers2Cells(
        Ma,
        nmark,
        PC_th,
        PC,
        weight_th,
        weight,
        x,
        y,
        Δ,
        param,
        param2;
        avgm = :arith,
    )

Interpolate a marker property to the extended cell-centered grid using
weighted bilinear interpolation.

Each active marker contributes to the four surrounding cell centers
according to its relative position within the corresponding grid cell.
Thread-local property and weight arrays are used to accumulate the marker
contributions in parallel. The thread-local arrays are subsequently summed
and normalized to obtain the interpolated cell-centered field.

The function can interpolate either marker temperature directly or a
phase-dependent material property specified through `param2`.

# Arguments

- `Ma`: Marker structure containing marker coordinates, phase IDs, and,
  where applicable, marker temperatures.
- `nmark`: Total number of markers.
- `PC_th`: Collection of thread-local property arrays defined on the
  extended cell-centered grid.
- `PC`: Extended cell-centered array receiving the interpolated property.
- `weight_th`: Collection of thread-local interpolation-weight arrays.
- `weight`: Extended cell-centered array receiving the accumulated weights.
- `x`: Horizontal coordinates of the extended cell centers.
- `y`: Vertical coordinates of the extended cell centers.
- `Δ`: Grid spacing in the horizontal and vertical directions.
- `param`: Property to interpolate:
  - `:thermal` interpolates `Ma.T`.
  - `:phase` interpolates the phase-dependent values contained in `param2`.
- `param2`: Collection containing the material-property value assigned to
  each phase. Phase ID `p` accesses the value `param2[p+1]`. This argument
  is not used when `param == :thermal`.

# Keyword arguments

- `avgm=:arith`: Averaging method used for phase-dependent properties:
  - `:arith`: Weighted arithmetic average.
  - `:harm`: Weighted harmonic average.
  - `:geom`: Weighted geometric average.

For thermal interpolation, the marker temperatures are combined using the
weighted arithmetic average.

# Notes

The arrays in `PC_th` and `weight_th` must have the same dimensions as
`PC` and `weight`, respectively, and their number must be compatible with
the available Julia threads.

Grid points receiving no marker contributions retain their previous values
from `PC`. The function modifies `PC`, `weight`, `PC_th`, and `weight_th`
in place.

# Example

```julia
Markers2Cells(
    Ma,
    nmark,
    MAVG.PC_th,
    D.ρ_ex,
    MAVG.wte_th,
    D.wte,
    x,
    y,
    Δ,
    :phase,
    ρ;
    avgm = :arith,
)

D.ρ .= D.ρ_ex[2:end-1, 2:end-1]
```
"""
@views function Markers2Cells(Ma,nmark,PC_th,PC,weight_th,weight,x,y,Δ,param,param2;avgm=:arith, verbose=false)
    PC0     =   copy(PC)
    PC      .*=     0.0
    weight  .*=     0.0
    if param==:thermal
        PM  =       copy(Ma.T)
        chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(PC_th[tid], 0)
                fill!(weight_th[tid], 0)
                for k in chunk
                    # Distance to the upper right corner ---
                    # Get the column:
                    dstx = Ma.x[k] - x.ce[1]
                    i = ceil(Int64, dstx / Δ.x) + 1                   
                    # Get the line:
                    dsty = Ma.y[k] - y.ce[1]
                    j = ceil(Int64, dsty /  Δ.y) + 1 
                    # Relative distances
                    Δxm = abs(x.ce[i] - Ma.x[k])/Δ.x
                    Δym = abs(y.ce[j] - Ma.y[k])/Δ.y
                    # Increment cell counts
                    PC_th[tid][i-1,j-1] += PM[k] * Δxm * Δym
                    PC_th[tid][i  ,j-1] += PM[k] * (1.0 - Δxm) * Δym
                    PC_th[tid][i-1,j  ] += PM[k] * Δxm * (1.0 - Δym)
                    PC_th[tid][i  ,j  ] += PM[k] * (1.0 - Δxm) * (1.0 - Δym)
                    
                    weight_th[tid][i-1,j-1]    += Δxm * Δym
                    weight_th[tid][i  ,j-1]    += (1.0 - Δxm) * Δym
                    weight_th[tid][i-1,j  ]    += Δxm * (1.0 - Δym)
                    weight_th[tid][i  ,j  ]    += (1.0 - Δxm) * (1.0 - Δym)
                end
            end
        end
    elseif param==:phase
        PM  =       copy(Ma.phase)
        chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(PC_th[tid], 0)
                fill!(weight_th[tid], 0)
                for k in chunk
                    # Distance to the upper right corner ---
                    # Get the column:
                    dstx = Ma.x[k] - x.ce[1]
                    i = ceil(Int64, dstx / Δ.x) + 1                   
                    # Get the line:
                    dsty = Ma.y[k] - y.ce[1]
                    j = ceil(Int64, dsty /  Δ.y) + 1 
                    # Relative distances
                    Δxm = abs(x.ce[i] - Ma.x[k])/Δ.x
                    Δym = abs(y.ce[j] - Ma.y[k])/Δ.y
                    # Increment cell counts
                    if avgm==:arith
                        PC_th[tid][i-1,j-1] += param2[PM[k]+1] * Δxm * Δym
                        PC_th[tid][i  ,j-1] += param2[PM[k]+1] * (1.0 - Δxm) * Δym
                        PC_th[tid][i-1,j  ] += param2[PM[k]+1] * Δxm * (1.0 - Δym)
                        PC_th[tid][i  ,j  ] += param2[PM[k]+1] * (1.0 - Δxm) * (1.0 - Δym)
                    elseif avgm==:harm
                        PC_th[tid][i-1,j-1] += (Δxm * Δym) / param2[PM[k]+1] 
                        PC_th[tid][i  ,j-1] += ((1.0 - Δxm) * Δym) / param2[PM[k]+1]
                        PC_th[tid][i-1,j  ] += (Δxm * (1.0 - Δym)) / param2[PM[k]+1]
                        PC_th[tid][i  ,j  ] += ((1.0 - Δxm) * (1.0 - Δym)) / param2[PM[k]+1]
                    elseif avgm==:geom
                        PC_th[tid][i-1,j-1] += log(param2[PM[k]+1]) * (Δxm * Δym)
                        PC_th[tid][i  ,j-1] += log(param2[PM[k]+1]) * ((1.0 - Δxm) * Δym)
                        PC_th[tid][i-1,j  ] += log(param2[PM[k]+1]) * (Δxm * (1.0 - Δym))
                        PC_th[tid][i  ,j  ] += log(param2[PM[k]+1]) * ((1.0 - Δxm) * (1.0 - Δym))
                    end
                    weight_th[tid][i-1,j-1]    += Δxm * Δym
                    weight_th[tid][i  ,j-1]    += (1.0 - Δxm) * Δym
                    weight_th[tid][i-1,j  ]    += Δxm * (1.0 - Δym)
                    weight_th[tid][i  ,j  ]    += (1.0 - Δxm) * (1.0 - Δym)
                end
            end
        end
    end
    
    PC      .= reduce(+, PC_th)
    weight  .= reduce(+, weight_th)
    
    if avgm==:arith
        PC      ./= weight
    elseif avgm ==:harm
        @. PC   =   weight / PC
    elseif avgm==:geom
        @. PC   =   exp(PC/weight) # PC^(1/weight)
    end

    if sum(isnan.(PC))>0
        if verbose 
            @printf("%i number(s) of cells without markers\n", sum(isnan.(PC)))
        end
        PC[isnan.(PC)]     .=  PC0[isnan.(PC)]
    end

    return
end

"""
    Markers2Vertices(
        Ma,
        nmark,
        PG_th,
        PG,
        weight_th,
        weight,
        x,
        y,
        Δ,
        param,
        param2;
        avgm = :arith,
    )

Interpolate a marker property to the grid vertices using weighted bilinear
interpolation.

Each marker contributes to the four surrounding vertices according to its
relative position within the corresponding grid cell. Thread-local property
and weight arrays are used to accumulate the marker contributions in
parallel. The thread-local arrays are subsequently summed and normalized to
obtain the interpolated vertex field.

The function can interpolate either marker temperature directly or a
phase-dependent material property specified through `param2`.

# Arguments

    - `Ma`: Marker structure containing marker coordinates, phase IDs, and,
    where applicable, marker temperatures.
    - `nmark`: Total number of markers.
    - `PG_th`: Collection of thread-local property arrays defined at the grid
    vertices.
    - `PG`: Vertex-centered array receiving the interpolated property.
    - `weight_th`: Collection of thread-local interpolation-weight arrays
    defined at the vertices.
    - `weight`: Vertex-centered array receiving the accumulated interpolation
    weights.
    - `x`: Horizontal coordinates of the grid vertices.
    - `y`: Vertical coordinates of the grid vertices.
    - `Δ`: Grid spacing in the horizontal and vertical directions.
    - `param`: Property to interpolate:
    - `:thermal` interpolates `Ma.T`.
    - `:phase` interpolates the phase-dependent values contained in `param2`.
    - `param2`: Collection containing the material-property value assigned to
    each phase. Phase ID `p` accesses the value `param2[p+1]`. This argument
    is not used when `param == :thermal`.

# Keyword arguments

    - `avgm=:arith`: Averaging method used for phase-dependent properties:
    - `:arith`: Weighted arithmetic average.
    - `:harm`: Weighted harmonic average.
    - `:geom`: Weighted geometric average.

For thermal interpolation, marker temperatures are combined using the
weighted arithmetic average.

# Notes

The arrays in `PG_th` and `weight_th` must have the same dimensions as
`PG` and `weight`, respectively, and their number must be compatible with
the available Julia threads.

Vertices receiving no marker contributions retain their previous values
from `PG`. The function modifies `PG`, `weight`, `PG_th`, and `weight_th`
in place.

# Example

```julia
Markers2Vertices(
    Ma,
    nmark,
    MAVG.PV_th,
    D.ηv,
    MAVG.wtv_th,
    D.wtv,
    x,
    y,
    Δ,
    :phase,
    η;
    avgm = :arith,
)
```
"""
@views function Markers2Vertices(Ma,nmark,PG_th,PG,weight_th,weight,x,y,Δ,param,param2;avgm=:arith, verbose=false)
    PG0     =   copy(PG)
    PG      .*=     0.0
    weight  .*=     0.0
    if param==:thermal
        PM  =       copy(Ma.T)
        chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(PG_th[tid], 0)
                fill!(weight_th[tid], 0)
                for k in chunk
                    # Get upper left corner ---
                    # Get the column:
                    dstx = Ma.x[k] - x.v[1]
                    i = ceil(Int64, dstx / Δ.x )
                    # Get the line:
                    dsty = Ma.y[k] - y.v[1]
                    j = ceil(Int64, dsty /  Δ.y ) + 1
                    # Relative distances
                    Δxm = abs(x.v[i] - Ma.x[k])/ Δ.x
                    Δym = abs(y.v[j] - Ma.y[k])/ Δ.y
                    # Increment cell counts
                    PG_th[tid][i  ,j  ]    += PM[k] * (1.0 - Δxm) * (1.0 - Δym)
                    PG_th[tid][i  ,j-1]    += PM[k] * (1.0 - Δxm) * Δym
                    PG_th[tid][i+1,j  ]    += PM[k] * Δxm * (1.0 - Δym)
                    PG_th[tid][i+1,j-1]    += PM[k] * Δxm * Δym
                    
                    weight_th[tid][i  ,j  ] += (1.0 - Δxm) * (1.0 - Δym)
                    weight_th[tid][i  ,j-1] += (1.0 - Δxm) * Δym
                    weight_th[tid][i+1,j  ] += Δxm * (1.0 - Δym)
                    weight_th[tid][i+1,j-1] += Δxm * Δym
                end
            end
        end
    elseif param==:phase
        PM  =       copy(Ma.phase)
        chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
        @sync for chunk in chunks
            @spawn begin
                tid = threadid()
                fill!(PG_th[tid], 0)
                fill!(weight_th[tid], 0)
                for k in chunk
                    # Get upper left corner ---
                    # Get the column:
                    dstx = Ma.x[k] - x.v[1]
                    i = ceil(Int64, dstx / Δ.x )
                    # Get the line:
                    dsty = Ma.y[k] - y.v[1]
                    j = ceil(Int64, dsty /  Δ.y ) + 1
                    # Relative distances
                    Δxm = abs(x.v[i] - Ma.x[k])/ Δ.x
                    Δym = abs(y.v[j] - Ma.y[k])/ Δ.y
                    # Increment cell counts
                    if avgm==:arith
                        PG_th[tid][i  ,j  ]    += param2[PM[k]+1] * (1.0 - Δxm) * (1.0 - Δym)
                        PG_th[tid][i  ,j-1]    += param2[PM[k]+1] * (1.0 - Δxm) * Δym
                        PG_th[tid][i+1,j  ]    += param2[PM[k]+1] * Δxm * (1.0 - Δym)
                        PG_th[tid][i+1,j-1]    += param2[PM[k]+1] * Δxm * Δym
                    elseif avgm==:harm
                        PG_th[tid][i  ,j  ]    += ((1.0 - Δxm) * (1.0 - Δym)) / param2[PM[k]+1]
                        PG_th[tid][i  ,j-1]    += ((1.0 - Δxm) * Δym) / param2[PM[k]+1]
                        PG_th[tid][i+1,j  ]    += (Δxm * (1.0 - Δym)) / param2[PM[k]+1]
                        PG_th[tid][i+1,j-1]    += (Δxm * Δym) / param2[PM[k]+1]
                    elseif avgm==:geom
                        PG_th[tid][i  ,j  ]    += log(param2[PM[k]+1]) * (1.0 - Δxm) * (1.0 - Δym)
                        PG_th[tid][i  ,j-1]    += log(param2[PM[k]+1]) * (1.0 - Δxm) * Δym
                        PG_th[tid][i+1,j  ]    += log(param2[PM[k]+1]) * Δxm * (1.0 - Δym)
                        PG_th[tid][i+1,j-1]    += log(param2[PM[k]+1]) * Δxm * Δym
                    end
                    weight_th[tid][i  ,j  ] += (1.0 - Δxm) * (1.0 - Δym)
                    weight_th[tid][i  ,j-1] += (1.0 - Δxm) * Δym
                    weight_th[tid][i+1,j  ] += Δxm * (1.0 - Δym)
                    weight_th[tid][i+1,j-1] += Δxm * Δym
                end
            end
        end
    end
    
    PG      .= reduce(+, PG_th)
    weight  .= reduce(+, weight_th)

    if avgm==:arith
        PG      ./= weight
    elseif avgm ==:harm
        @. PG   =   weight / PG
    elseif avgm==:geom
        @. PG   =   exp(PG/weight) # PC^(1/weight)
    end

    if sum(isnan.(PG))>0
        if verbose 
            @printf("%i number(s) of vertices without markers\n",sum(isnan.(PG)))
        end
        PG[isnan.(PG)]     .=  PG0[isnan.(PG)]
    end

    return
end

"""
    AdvectTracer2D(
        Ma,
        nmark,
        D,
        x,
        y,
        dt,
        Δ,
        NC,
        rkw,
        rkv;
        style = 1,
    )

Advect two-dimensional Lagrangian tracers using a fourth-order Runge–Kutta
scheme.

For each active tracer, the velocity is interpolated from the Eulerian grid
at the intermediate Runge–Kutta positions. The stage velocities are combined
using the coefficients in `rkw`, while `rkv` defines the intermediate
positions used during the four Runge–Kutta stages.

Only tracers with a non-negative phase ID are advected.

# Arguments

- `Ma`: Marker structure containing the tracer coordinates and phase IDs.
- `nmark`: Total number of tracers.
- `D`: Data structure containing the staggered and cell-centered velocity
  fields.
- `x`: Horizontal coordinates of the staggered and cell-centered grids.
- `y`: Vertical coordinates of the staggered and cell-centered grids.
- `dt`: Time-step size.
- `Δ`: Grid spacing in the horizontal and vertical directions.
- `NC`: Number of cell centers in each coordinate direction.
- `rkw`: Runge–Kutta weights used to combine the stage velocities.
- `rkv`: Runge–Kutta coefficients used to calculate the intermediate tracer
  positions.

# Keyword arguments

- `style=1`: Velocity-interpolation method:
  - `1`: Bilinear interpolation from the staggered `vₓ` and `vᵧ` nodes.
  - `2`: Combined interpolation using staggered and cell-centered
    velocities.
  - `3`: Interpolation from the staggered velocity nodes with the optional
    higher-order correction enabled.

# Returns

Returns the modified marker structure `Ma`. The tracer coordinates `Ma.x`
and `Ma.y` are updated in place.

# Notes

Markers with

```julia
Ma.phase[k] < 0
```

are treated as inactive and are not advected.

The velocity-interpolation routines are evaluated at the intermediate
Runge–Kutta positions by temporarily updating the tracer coordinates during
each stage.

# Example

rkw = (1.0 / 6.0) .* [1.0, 2.0, 2.0, 1.0]
rkv = (1.0 / 2.0) .* [1.0, 1.0, 2.0, 2.0]

AdvectTracer2D(
    Ma,
    nmark,
    D,
    x,
    y,
    dt,
    Δ,
    NC,
    rkw,
    rkv;
    style = 1,
)
```
"""
@views function AdvectTracer2D(Ma,nmark,D,x,y,dt,Δ,NC,rkw,rkv;style=1)
    @threads for k = 1:nmark
        if (Ma.phase[k]>=0)
            x0  =   Ma.x[k]
            y0  =   Ma.y[k]
            vx  =   0.0
            vy  =   0.0
            # Runge-Kutta loop ---
            for rk=1:4
                # Interp velocity from grid ---
                if style == 1 # Bilinear velocity interp (original is Markers_divergence_ALLSCHEMES_RK4.m)
                    vxm = VxFromVxNodes(D.vx, k, Ma, x, y, Δ, NC, 0)
                    vym = VyFromVyNodes(D.vy, k, Ma, x, y, Δ, NC, 0)
                elseif style == 2
                    vxx = VxFromVxNodes(D.vx, k, Ma, x, y, Δ, NC, 0)
                    vyy = VyFromVxNodes(D.vy, k, Ma, x, y, Δ, NC, 0)
                    vxp, vyp = VxVyFromPrNodes(D.vxc ,D.vyc, k, Ma, x, y, Δ, NC)
                    vxm = itpw*vxp + (1.0-itpw)*vxx
                    vym = itpw*vyp + (1.0-itpw)*vyy
                elseif style == 3
                    vxm = VxFromVxNodes(D.vx, k, Ma, x, y, Δ, NC, 1)
                    vym = VyFromVxNodes(D.vy, k, Ma, x, y, Δ, NC, 1)
                end
                # Temporary RK advection steps ---
                Ma.x[k]     =   x0 + rkv[rk]*dt*vxm
                Ma.y[k]     =   y0 + rkv[rk]*dt*vym
                # Average final velocity ---
                vx    += rkw[rk]*vxm
                vy    += rkw[rk]*vym
            end
            # Advect points ---
            Ma.x[k]     =   x0 + rkv[4]*dt*vx
            Ma.y[k]     =   y0 + rkv[4]*dt*vy
        end
    end
    return Ma
end
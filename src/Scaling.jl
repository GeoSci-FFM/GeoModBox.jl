"""
    Constants()

Structure containing the characteristic scaling quantities used for
non-dimensionalization and post-processing.

The structure stores the reference scales for length, time, velocity,
stress, temperature, and volumetric heat production. By default, all
values are initialized to zero and can be assigned by the user depending
on the chosen scaling.

# Fields

    - `hsc`: Characteristic length scale.
    - `tsc`: Characteristic time scale.
    - `vsc`: Characteristic velocity scale.
    - `τsc`: Characteristic stress scale.
    - `Tsc`: Characteristic temperature scale.
    - `Qsc`: Characteristic volumetric heat-production scale.

# Default values

```julia
hsc = 0.0
tsc = 0.0
vsc = 0.0
τsc = 0.0
Tsc = 0.0
Qsc = 0.0
```
# Example

S = Constants()

S.hsc = 100e3      # m
S.tsc = 1.0e6      # yr
S.vsc = 0.1        # cm/yr
S.τsc = 1.0e8      # Pa
S.Tsc = 1300.0     # °C
S.Qsc = 1.0e-6     # W/m³
"""
@kwdef mutable struct Constants
    hsc ::Float64 = 0.0
    tsc ::Float64 = 0.0
    vsc ::Float64 = 0.0
    τsc ::Float64 = 0.0
    Tsc ::Float64 = 0.0
    Qsc ::Float64 = 0.0
end

"""
    ScalingConstants!(M, P)

Compute the characteristic scaling constants used for
non-dimensionalization.

The scaling is based on the model height and the thermal diffusion time,
which are commonly used as characteristic scales in mantle convection
problems. The resulting structure contains the reference scales for
length, time, velocity, stress, temperature, and volumetric heat
production.

The scaling constants are defined as

```math
h_{\\mathrm{sc}} = H,

t_{\\mathrm{sc}} = \\frac{H^2}{\\kappa},

v_{\\mathrm{sc}} = \\frac{\\kappa}{H},

\tau_{\\mathrm{sc}} = \\frac{\\eta_0 \\kappa}{H^2},

T_{\\mathrm{sc}} = \\Delta T,

Q_{\\mathrm{sc}} = \\frac{\\Delta T \\kappa \\rho_0 c_p}{H^2},
```

where
    - H = ymax - ymin is the model height,
    - κ is the thermal diffusivity,
    - η₀ is the reference viscosity,
    - ΔT is the reference temperature difference,
    - ρ₀ is the reference density, and
    - cₚ is the specific heat capacity.

# Arguments

    - M: Model geometry.
    - P: Structure containing the physical parameters.

# Returns

Returns a Constants structure containing the characteristic scaling constants in SI units.

# Example

M = Geometry(0, 1000, -1000, 0)
P = Physics()
S = ScalingConstants!(M, P)

S.vsc

# Notes

The returned scaling constants can be used to convert dimensional model parameters 
to their non-dimensional counterparts, or to convert non-dimensional simulation 
results back to SI units for analysis and visualization. 
""" 
function ScalingConstants!(M,P)

    S       =   Constants()

    S.hsc   =   (M.ymax-M.ymin)                         #   Length scale [ m ]
    S.tsc   =   (M.ymax-M.ymin)^2 / P.κ                 #   Time scale [ s ]
    S.vsc   =   P.κ / (M.ymax-M.ymin)                   #   Velocity scale [ m/s ]
    S.τsc   =   (P.η₀ * P.κ)/(M.ymax-M.ymin)^2          #   Stress scale [ Pa ]
    S.Tsc   =   P.ΔT                                    #   Temperature scale [ K ]
    S.Qsc   =   (P.ΔT*P.κ*P.ρ₀*P.cp)/(M.ymax-M.ymin)^2  #   Heat source scale [ w/m³ ]

    return S
    
end

"""
    ScaleParameters!(S, M, Δ, T, P, D)

Scale the model geometry, time parameters, physical parameters, and field
variables using the characteristic scales stored in `S`.

The function converts the dimensional model setup to its non-dimensional
form. The characteristic scales are typically obtained from
`ScalingConstants!()` and are based on the model height and the thermal
diffusion time. The scaling is performed in place by modifying the
corresponding structures.

The following quantities are non-dimensionalized:

### Geometry

```math
x=\\frac{x}{h_{\\mathrm{sc}}}, \\qquad
y=\\frac{y}{h_{\\mathrm{sc}}}
```

including the model boundaries and the grid spacing.

### Time

```math
t=\\frac{t}{t_{\\mathrm{sc}}}
```

including the maximum simulation time and all time-step estimates.

### Temperature

The boundary temperatures are converted from Kelvin to non-dimensional
temperature using

```math
T=\\frac{T-273.15}{T_{\\mathrm{sc}}}.
```

### Volumetric heat production

```math
Q=\\frac{Q}{Q_{\\mathrm{sc}}}.
```

The reference viscosity scaling

```math
\\eta=\\frac{\\eta}{\\eta_0}
```

can optionally be applied if viscosity is stored in dimensional units.

# Arguments

    - `S`: Structure containing the characteristic scaling constants.
    - `M`: Model geometry.
    - `Δ`: Grid spacing.
    - `T`: Time-parameter structure.
    - `P`: Physical-parameter structure.
    - `D`: Structure containing the model fields.

# Notes

The function modifies all supplied structures in place. After calling
`ScaleParameters!`, the governing equations can be solved in their
non-dimensional form.

Typically, the scaling workflow is

```julia
S = ScalingConstants!(M, P)
ScaleParameters!(S, M, Δ, T, P, D)
```

# Example

```julia
M = Geometry(0, 1000, -1000, 0)
P = Physics()

S = ScalingConstants!(M, P)

ScaleParameters!(S, M, Δ, T, P, D)
```
"""
function ScaleParameters!(S,M,Δ,T,P,D)
    # Geometry ---
    M.xmin   /=  S.hsc
    M.xmax   /=  S.hsc
    M.ymin   /=  S.hsc
    M.ymax   /=  S.hsc
    Δ.x      /=  S.hsc
    Δ.y      /=  S.hsc         
    # Time --- 
    T.tmax   /=  S.tsc
    T.Δc     /=  S.tsc
    T.Δd     /=  S.tsc
    T.Δ      /=  S.tsc
    # Temperature etc. ---
    P.Ttop   =   (P.Ttop - 273.15)/S.Tsc
    P.Tbot   =   (P.Tbot - 273.15)/S.Tsc
    @. D.Q      /=  S.Qsc
    # # Viscosity ---
    # @. D.η      /=  P.η₀
end
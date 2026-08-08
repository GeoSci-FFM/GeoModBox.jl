# Mixed Heated Convection

This directory contains examples of several types of mixed heated thermal convection, each driven by distinct heating mechanisms. 

# Isoviscous Thermal Convection

- [Purely internally heated](InternallyHeated.jl)
- [Purely basally heated](BottomHeated.jl)
- [Internally and basally heated](MixedHeated.jl)

---

# Temperature Dependent Viscosity 

- [Purely basally heated](BottomHeated_VarEta.jl)
- [Blankenbach benchmark](Blankenbach_var_eta.jl)

--- 

All simulations use non-dimensional quantities. The governing equations are solved as follows:
- The momentum conservation equation is solved using the general solver.
- The energy conservation equation is solved using the general solver with a Crank–Nicolson discretization.
- Temperature advection is performed using a semi-Lagrangian scheme.

For further details on the numerical methods used for each governing equation and for additional information on each example, please refer to the [GeoModBox.jl documentation](https://geosci-ffm.github.io/GeoModBox.jl/).

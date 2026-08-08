# Thermo-mechanical Shear Localization

This directory contains an example for the thermo-mechanical shear localization in a viscous medium due to viscous dissipation (shear heating). The files included are: 

- A script to generate the data of a [model series for the viscous including](GenerateData.jl) model following the setup of Duretz et al. (2014)

- Script which [plots the data](PlotData.jl) generated in GenerateData.jl

- [Main script](ShearHeatingShearBands.jl) to solve the thermo-mechanical shear localization model follwing the setup of Duretz et al. (2014)

For further details on the discretization methods and implementation of each example, please refer to the [GeoModBox.jl documentation](https://geosci-ffm.github.io/GeoModBox.jl/).

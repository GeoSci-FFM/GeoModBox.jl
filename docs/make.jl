using Documenter
using GeoModBox

GMB_root_dir = dirname(@__DIR__)

license = read(joinpath(GMB_root_dir, "LICENSE"), String)
write(joinpath(@__DIR__, "src", "man", "license.md"), license)

makedocs(
    sitename = "GeoModBox.jl",
    format = Documenter.HTML(),
    warnonly = [:missing_docs],
    modules = [GeoModBox],
    pages = [
        "Home" => "index.md", 
        "Installation" => "man/Installation.md",
        "Governing Equation" => Any[
            "Solution" => Any[
                "General" => "man/theory/GESolution.md", 
                "Specifications" => "man/Ini.md", 
            ],
            "Heat Diffusion Equation" => Any[
                "General" => "man/theory/DiffMain.md",
                "1D" => "man/theory/DiffOneD.md",
                "2D" => "man/theory/DiffTwoD.md",
            ],
            "Advection Equation" => Any[
                "General" => "man/theory/AdvectMain.md",
                "1D" => "man/theory/AdvOneD.md",
                "2D" => "man/theory/AdvTwoD.md",
            ],
            "Momentum Equation" => Any[
                "General" => "man/theory/MomentumMain.md",
                "1D" => "man/theory/MomentumOneD.md",
                "2D" => "man/theory/MomentumTwoD.md",
            ],
        ],
        "Exercises" => Any[
            "General" => "man/exercises/Exercises.md",
            "01 - Euler Advection" => "man/exercises/01_Euler_Advection.md",
            "02 - 1D Heat Diffusion (explicit)" => "man/exercises/02_1D_Heat_explicit.md",
            "03 - 1D Heat Diffusion (implicit)" => "man/exercises/03_1D_Heat_implicit.md",
            "04 - 2D Heat Diffusion (stationary)" => "man/exercises/04_2D_Diffusion_Stationary.md",
            "05a - 2D Heat Diffusion (Plume)" => "man/exercises/05_2D_Diffusion_TD_Plume.md", 
            "05b - 2D Heat Diffusion (Sill)" => "man/exercises/05_2D_Diffusion_TD_Sill.md", 
            "06 - 1D Advection" => "man/exercises/06_1D_Advection.md", 
            "07 - 2D Energy Conservation" => "man/exercises/07_2D_Energy_Equation.md",
            "08 - 1D Stokes" => "man/exercises/08_1D_Stokes.md",
            "09 - 2D Falling Block (steady state)" => "man/exercises/09_2D_Falling_Block.md", 
            "10 - 2D Falling Block (time-dep)" => "man/exercises/10_2D_Falling_Block_td.md", 
            "11 - 2D Thermal Convection" => "man/exercises/11_2D_Thermal_Convection.md", 
            "12 - 2D Thermal Convection (scaled)" => "man/exercises/12_2D_Thermal_Convection_scaled.md", 
            "13 - Blankenbach Benchmark" => "man/exercises/13_Blankenbach_Benchmark.md", 
        ],
        "Examples" => Any[
            "General" => "man/examples/Examples.md",
            "Diffusion Equation" => Any[
                "Oceanic Geotherm (1D)" => "man/examples/Diffusion/OceanicGeotherm.md", 
                "Continental Geotherm (1D)" => "man/examples/Diffusion/ContinentalGeotherm.md", 
                "Gaussian Diffusion (1D)" => "man/examples/Diffusion/GaussianDiffusion1D.md", 
                "Backward Euler (2D)" => "man/examples/Diffusion/BackwardEuler_DC.md", 
                "Forward Euler (2D)" => "man/examples/Diffusion/ForwardEuler_DC.md", 
                "Gaussian Diffusion (2D)" => "man/examples/Diffusion/GaussianDiffusion2D.md", 
                "Poisson Problem (2D)" => "man/examples/Diffusion/PoissonRestest.md", 
                "Poisson Problem; variable k (2D)" => "man/examples/Diffusion/PoissonVariablek.md", 
            ],
            "Advection Equtaion" => Any[
                "Advection (2D)" => "man/examples/Advection/Advection2D.md",
                "Advection Resolution Test (2D)" => "man/examples/Advection/AdvectionRestest2D.md",
            ],
            "Stokes Equation" => Any[
                "Channel Flow (1D)" => "man/examples/Stokes/ChannelFlow1D.md", 
                "Falling Block Benchmark (Gerya, 2019)" => "man/examples/Stokes/FallingBlockBenchmark.md", 
                "Falling Block" => "man/examples/Stokes/FallingBlockDC.md", 
                "Rayleigh Taylor Instability (RTI)" => "man/examples/Stokes/RTI.md", 
                "RTI - Growth Rate (Ramberg, 1968)" => "man/examples/Stokes/RTI_growth_rate.md", 
                "RTI (Van Keken et al., 1997)" => "man/examples/Stokes/VanKekenBenchmark.md",
                "Viscous Inclusion (Schmid, 2002)" => "man/examples/Stokes/ViscousInclusion.md",
            ],
            "Mixed Thermal Convection" => Any[
                "Overview" => "man/examples/Convection/Overview_Convection.md",
                "Bottom Heated" => "man/examples/Convection/BottomHeatedConvection.md",
                "Internally Heated" => "man/examples/Convection/InternallyHeatedConvection.md",
                "Mixed Heated" => "man/examples/Convection/MixedHeatedConvection.md",
                "Bottom Heated - Variable Visc" => "man/examples/Convection/BottomHeatedConvection_VE.md",
                "Blankenbach et al. (1989) - Variable Visc" => "man/examples/Convection/Blankenbach_Var_Eta.md",
            ],
            "Thermo-Mechanical Shear Localization" => Any[
                "Viscous shear localization (Duretz et al., 2014)" => "man/examples/StrainLocalization/ShearBands.md",
            ]
        ],
        "List of functions" => "man/listoffunctions.md",
        "License" => "man/license.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
withenv("GITHUB_REPOSITORY" => "GeoSci-FFM/GeoModBox.jl") do
    deploydocs(
        # repo = "https://lukasfuchs.github.io/GeoModBox.jl",
        repo = "github.com/GeoSci-FFM/GeoModBox.jl",
        devbranch = "main",
        push_preview = true
    )
end
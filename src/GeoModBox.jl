module GeoModBox

    using Statistics: mean
    using ExtendableSparse, LinearAlgebra

    module HeatEquation

        module TwoD
            # Handle analytical solution for 2D Diffusion2D_Gaussian ---
            include("./HeatEquation/AnalyticsDiffusion2D.jl")
            export AnalyticalSolution2D!, BoundaryConditions2D!	    

            # Implicit solver ---
            include("./HeatEquation/BackwardEuler.jl")
            export ComputeResiduals2D!, AssembleMatrix2D, 
                    BackwardEuler2Dc!

            # Poisson solver ---
            include("./HeatEquation/PoissonSolvers.jl")
            export Poisson2Dc!, Poisson2D!

            # Explicit solver ---
            include("./HeatEquation/ForwardEuler.jl")
            export ForwardEuler2Dc!

            # Cranck-Nicolson Approach ---
            include("./HeatEquation/CNA.jl")
            export CNA2Dc! 

            # Alternate Direct Implicit Method ---
            include("./HeatEquation/ADI.jl")
            export ADI2Dc!

        end

        module OneD
            # 1D solver ---
            include("./HeatEquation/1Dsolvers.jl")
            export ForwardEuler1Dc!, BackwardEuler1Dc!,
                ComputeResiduals1Dc!, AssembleMatrix1Dc!, CNA1Dc!,
                ForwardEuler1D!
        end                
    end

    module AdvectionEquation
    
        module OneD
            # 1D solver ---
            include("./AdvectionEquation/1Dsolvers.jl")
            export upwind1D!, lax1D!, slf1D!, semilag1D!, RK4O1D!

            # Tracer options ---
            include("./Tracers/ItpTracers.jl")
            export Itp1D_Centers2Markers!, Itp1D_Markers2Centers!
        end

        module TwoD
            # 2D solver --- 
            include("./AdvectionEquation/2Dsolver.jl")
            export upwindc2D!, slfc2D!, semilagc2D!
        end
    end

    module InitialCondition
        # Initial Conditions ---
        include("./InitialCondition/2Dini.jl")        
        export IniVelocity!, IniTemperature!
    end
end 

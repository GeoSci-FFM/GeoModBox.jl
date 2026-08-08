using GeoModBox.MomentumEquation.OneD
using ExtendableSparse

@testset "1-D Stokes equation" begin

    # ================================================================ #
    # Couette flow
    #
    # Constant viscosity, zero pressure gradient, and prescribed
    # velocities at the lower and upper boundaries:
    #
    #     vx(0) = 0
    #     vx(H) = V
    #
    # The analytical solution is
    #
    #     vx(y) = V * y / H
    #
    # ================================================================ #

    @testset "Couette flow - direct solver" begin

        nc = 16
        H  = 1.0
        Δy = H / nc

        η0 = 2.0
        V  = 1.0

        # Velocity is defined at cell centroids
        yc = range(Δy / 2, H - Δy / 2; length=nc)

        # Viscosity is defined at cell faces
        η = fill(η0, nc + 1)

        BC = (
            type = (
                S = :Dirichlet,
                N = :Dirichlet,
            ),
            val = (
                S = 0.0,
                N = V,
            ),
        )

        # Zero pressure gradient
        rhs = zeros(nc)

        vx = zeros(nc)

        K = ExtendableSparseMatrix(nc, nc)

        Stokes_1D_direct(
            vx,
            η,
            Δy,
            nc,
            BC,
            K,
            rhs,
        )

        # Analytical Couette-flow solution
        vexact = @. V * yc / H

        @test vx ≈ vexact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Residual of analytical Couette flow
    # ================================================================ #

    @testset "Couette flow - residual" begin

        nc = 16
        H  = 1.0
        Δy = H / nc

        η0 = 2.0
        V  = 1.0

        yc = range(Δy / 2, H - Δy / 2; length=nc)

        BC = (
            type = (
                S = :Dirichlet,
                N = :Dirichlet,
            ),
            val = (
                S = 0.0,
                N = V,
            ),
        )

        # Analytical velocity at the cell centroids
        vx = @. V * yc / H

        # Minimal data structure required by ComputeStokesResiduals1D!
        D = (
            vₓ      = copy(vx),
            vₓₑ     = zeros(nc + 2),
            η       = fill(η0, nc + 1),
            ∂vₓ∂y   = zeros(nc + 1),
            τxy     = zeros(nc + 1),
            ∂τxy∂y  = zeros(nc),
            R       = zeros(nc),
        )

        ∂P∂x = 0.0

        ComputeStokesResiduals1D!(
            D,
            ∂P∂x,
            Δy,
            BC,
        )

        # Linear velocity -> constant strain rate
        @test D.∂vₓ∂y ≈ fill(V / H, nc + 1) atol=1e-12 rtol=1e-12

        # Constant viscosity -> constant shear stress
        @test D.τxy ≈ fill(η0 * V / H, nc + 1) atol=1e-12 rtol=1e-12

        # Constant shear stress -> zero stress divergence
        @test D.∂τxy∂y ≈ zeros(nc) atol=1e-12 rtol=1e-12

        # Therefore the Stokes residual must vanish
        @test D.R ≈ zeros(nc) atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Singular boundary-condition combination
    # ================================================================ #

    @testset "Reject double-Neumann system" begin

        nc = 8
        Δy = 1.0 / nc

        η   = ones(nc + 1)
        vx  = zeros(nc)
        rhs = zeros(nc)

        BC = (
            type = (
                S = :Neumann,
                N = :Neumann,
            ),
            val = (
                S = 0.0,
                N = 0.0,
            ),
        )

        K = ExtendableSparseMatrix(nc, nc)

        @test_throws ArgumentError Stokes_1D_direct(
            vx,
            η,
            Δy,
            nc,
            BC,
            K,
            rhs,
        )

        K = ExtendableSparseMatrix(nc, nc)

        @test_throws ArgumentError AssembleStokesMatrix1D(
            nc,
            η,
            Δy,
            BC,
            K,
        )

    end

end

@testset "Poiseuille flow" begin

    nc = 16
    H  = 1.0
    Δy = H / nc

    η0   = 2.0
    ∂P∂x = -4.0

    yc = range(Δy / 2, H - Δy / 2; length=nc)

    η = fill(η0, nc + 1)

    BC = (
        type = (
            S = :Dirichlet,
            N = :Dirichlet,
        ),
        val = (
            S = 0.0,
            N = 0.0,
        ),
    )

    rhs = fill(∂P∂x, nc)

    vx = zeros(nc)
    K  = ExtendableSparseMatrix(nc, nc)

    Stokes_1D_direct(
        vx,
        η,
        Δy,
        nc,
        BC,
        K,
        rhs,
    )

    # ------------------------------------------------------------ #
    # Continuous analytical solution
    # ------------------------------------------------------------ #

    vcont = @. (∂P∂x / (2 * η0)) * yc * (yc - H)

    # ------------------------------------------------------------ #
    # Exact solution of the cell-centered discrete system
    #
    # The Δy²/4 correction results from the ghost-node
    # implementation of the Dirichlet boundary conditions.
    # ------------------------------------------------------------ #

    vdisc = @. (∂P∂x / (2 * η0)) *
                 (yc * (yc - H) - Δy^2 / 4)

    @test vx ≈ vdisc atol=1e-12 rtol=1e-12


    # ================================================================ #
    # Residual of discrete solution
    # ================================================================ #

    D = (
        vₓ      = copy(vdisc),
        vₓₑ     = zeros(nc + 2),
        η       = copy(η),
        ∂vₓ∂y   = zeros(nc + 1),
        τxy     = zeros(nc + 1),
        ∂τxy∂y  = zeros(nc),
        R       = zeros(nc),
    )

    ComputeStokesResiduals1D!(
        D,
        ∂P∂x,
        Δy,
        BC,
    )

    @test D.R ≈ zeros(nc) atol=1e-12 rtol=1e-12


    # ================================================================ #
    # Continuous solution should converge with second-order accuracy
    # ================================================================ #

    error = maximum(abs.(vx .- vcont))

    @test error ≈ abs(∂P∂x / (2 * η0)) * Δy^2 / 4
end
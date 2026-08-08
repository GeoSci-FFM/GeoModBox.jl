using GeoModBox.AdvectionEquation.OneD

@testset "1-D advection" begin

    # ================================================================ #
    # Constant field
    #
    # A spatially constant scalar field has
    #
    #     ∂A/∂x = 0
    #
    # and must therefore remain unchanged during advection.
    # ================================================================ #

    @testset "Constant scalar field" begin

        nx  = 32
        Δx  = 1.0 / (nx - 1)
        vx  = 1.0
        CFL = 0.5
        Δt  = CFL * Δx / abs(vx)

        A0 = 5.0
        Aexact = fill(A0, nx)

        # ------------------------------------------------------------ #
        # Upwind
        # ------------------------------------------------------------ #

        A = copy(Aexact)

        upwind1D!(A, vx, Δt, Δx)

        @test A ≈ Aexact atol=1e-12 rtol=1e-12

        # Also test negative velocity because upwind1D! contains
        # separate branches for vx > 0 and vx < 0.

        A = copy(Aexact)

        upwind1D!(A, -vx, Δt, Δx)

        @test A ≈ Aexact atol=1e-12 rtol=1e-12

        # ------------------------------------------------------------ #
        # Lax-Friedrichs
        # ------------------------------------------------------------ #

        A = copy(Aexact)

        lax1D!(A, vx, Δt, Δx)

        @test A ≈ Aexact atol=1e-12 rtol=1e-12

        # ------------------------------------------------------------ #
        # Staggered leapfrog
        # ------------------------------------------------------------ #

        A     = copy(Aexact)
        Aold2 = copy(Aexact)

        slf1D!(A, Aold2, vx, Δt, Δx)

        @test A ≈ Aexact atol=1e-12 rtol=1e-12
        @test Aold2 ≈ Aexact atol=1e-12 rtol=1e-12

        # ------------------------------------------------------------ #
        # Semi-Lagrangian
        # ------------------------------------------------------------ #

        xc = range(0.0, 1.0; length=nx)
        A  = copy(Aexact)

        semilag1D!(A, xc, vx, Δt)

        @test A ≈ Aexact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # RK4 tracer advection
    #
    # For constant velocity, RK4 reduces to the exact translation
    #
    #     x_new = x_old + vx * Δt
    #
    # with periodic wrapping at the domain boundaries.
    # ================================================================ #

    @testset "RK4 tracer translation" begin

        xmin = 0.0
        xmax = 1.0

        vx = 0.2
        Δt = 0.5

        x = [0.10, 0.35, 0.70]

        xexact = x .+ vx * Δt

        RK4O1D!(x, Δt, vx, xmin, xmax)

        @test x ≈ xexact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # RK4 periodic wrapping
    # ================================================================ #

    @testset "RK4 periodic boundaries" begin

        xmin = 0.0
        xmax = 1.0

        vx = 0.4
        Δt = 0.5

        x = [0.10, 0.50, 0.90]

        # Translation by +0.2:
        #
        # 0.10 -> 0.30
        # 0.50 -> 0.70
        # 0.90 -> 1.10 -> 0.10 (periodic)
        xexact = [0.30, 0.70, 0.10]

        RK4O1D!(x, Δt, vx, xmin, xmax)

        @test x ≈ xexact atol=1e-12 rtol=1e-12

    end

end

@testset "Upwind exact one-cell translation" begin

    nx = 8
    Δx = 1.0
    vx = 1.0

    # CFL = 1
    Δt = Δx / abs(vx)

    # Boundary values are not updated by upwind1D!, so choose a field
    # for which we can test the interior nodes exactly.
    A0 = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]

    # ------------------------------------------------------------ #
    # Positive velocity
    # ------------------------------------------------------------ #

    A = copy(A0)

    upwind1D!(A, vx, Δt, Δx)

    # For vx > 0 and CFL = 1:
    #
    #     Aᵢⁿ⁺¹ = Aᵢ₋₁ⁿ
    #
    @test A[2:end-1] ≈ A0[1:end-2] atol=1e-12 rtol=1e-12

    # ------------------------------------------------------------ #
    # Negative velocity
    # ------------------------------------------------------------ #

    A = copy(A0)

    upwind1D!(A, -vx, Δt, Δx)

    # For vx < 0 and |CFL| = 1:
    #
    #     Aᵢⁿ⁺¹ = Aᵢ₊₁ⁿ
    #
    @test A[2:end-1] ≈ A0[3:end] atol=1e-12 rtol=1e-12

end
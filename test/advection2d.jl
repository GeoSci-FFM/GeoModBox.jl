using GeoModBox.AdvectionEquation.TwoD

@testset "2-D advection" begin

    # ================================================================ #
    # Constant scalar field
    #
    # A spatially constant property satisfies
    #
    #     ∂P/∂x = ∂P/∂y = 0
    #
    # and must therefore remain unchanged during advection.
    # ================================================================ #

    @testset "Constant scalar field" begin

        NC = (
            x = 8,
            y = 6,
        )

        Lx = 1.0
        Ly = 1.0

        Δx = Lx / NC.x
        Δy = Ly / NC.y

        Δt = 0.25 * min(Δx, Δy)

        P0 = 5.0

        Pexact = fill(P0, NC.x, NC.y)

        # Non-zero velocity in both directions
        vxc = fill(0.5, NC.x, NC.y)
        vyc = fill(-0.25, NC.x, NC.y)


        # ============================================================ #
        # Upwind
        # ============================================================ #

        P    = copy(Pexact)
        P_ex = fill(P0, NC.x + 2, NC.y + 2)

        upwindc2D!(
            P,
            P_ex,
            vxc,
            vyc,
            NC,
            Δt,
            Δx,
            Δy,
        )

        @test P ≈ Pexact atol=1e-12 rtol=1e-12


        # ============================================================ #
        # Staggered leapfrog
        # ============================================================ #

        P     = copy(Pexact)
        P_ex  = fill(P0, NC.x + 2, NC.y + 2)
        P_exo = fill(P0, NC.x + 2, NC.y + 2)

        slfc2D!(
            P,
            P_ex,
            P_exo,
            vxc,
            vyc,
            NC,
            Δt,
            Δx,
            Δy,
        )

        @test P ≈ Pexact atol=1e-12 rtol=1e-12


        # ============================================================ #
        # Semi-Lagrangian
        # ============================================================ #

        xc = range(Δx / 2, Lx - Δx / 2; length=NC.x)
        yc = range(Δy / 2, Ly - Δy / 2; length=NC.y)

        xce = range(-Δx / 2, Lx + Δx / 2; length=NC.x + 2)
        yce = range(-Δy / 2, Ly + Δy / 2; length=NC.y + 2)

        x2d = repeat(collect(xc), 1, NC.y)
        y2d = repeat(collect(yc)', NC.x, 1)

        x = (
            c   = xc,
            ce  = xce,
            c2d = x2d,
        )

        y = (
            c   = yc,
            ce  = yce,
            c2d = y2d,
        )

        P    = copy(Pexact)
        P_ex = fill(P0, NC.x + 2, NC.y + 2)

        # Previous velocity equal to current velocity
        vxo = copy(vxc)
        vyo = copy(vyc)

        semilagc2D!(
            P,
            P_ex,
            vxc,
            vyc,
            vxo,
            vyo,
            x,
            y,
            Δt,
        )

        @test P ≈ Pexact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Upwind exact one-cell translation
    #
    # For CFL = 1 in one coordinate direction, first-order upwind
    # reduces to an exact one-cell shift of the discrete field.
    # ================================================================ #

    @testset "Upwind exact one-cell translation" begin

        NC = (
            x = 8,
            y = 6,
        )

        Δx = 1.0
        Δy = 1.0

        # ------------------------------------------------------------ #
        # Translation in +x
        # ------------------------------------------------------------ #

        vx = 1.0
        vy = 0.0

        Δt = Δx / abs(vx)

        vxc = fill(vx, NC.x, NC.y)
        vyc = fill(vy, NC.x, NC.y)

        P0 = reshape(
            collect(1.0:(NC.x * NC.y)),
            NC.x,
            NC.y,
        )

        P_ex = zeros(NC.x + 2, NC.y + 2)
        P_ex[2:end-1,2:end-1] .= P0

        # Use constant extrapolation for the ghost nodes.
        P_ex[1,2:end-1]   .= P0[1,:]
        P_ex[end,2:end-1] .= P0[end,:]

        P_ex[2:end-1,1]   .= P0[:,1]
        P_ex[2:end-1,end] .= P0[:,end]

        P = copy(P0)

        upwindc2D!(
            P,
            P_ex,
            vxc,
            vyc,
            NC,
            Δt,
            Δx,
            Δy,
        )

        # For positive vx and CFLx = 1:
        #
        #     Pᵢⱼⁿ⁺¹ = Pᵢ₋₁,ⱼⁿ
        #
        @test P[2:end,:] ≈ P0[1:end-1,:] atol=1e-12 rtol=1e-12


        # ------------------------------------------------------------ #
        # Translation in -y
        # ------------------------------------------------------------ #

        vx = 0.0
        vy = -1.0

        Δt = Δy / abs(vy)

        vxc .= vx
        vyc .= vy

        P_ex .= 0.0
        P_ex[2:end-1,2:end-1] .= P0

        P_ex[1,2:end-1]   .= P0[1,:]
        P_ex[end,2:end-1] .= P0[end,:]

        P_ex[2:end-1,1]   .= P0[:,1]
        P_ex[2:end-1,end] .= P0[:,end]

        P = copy(P0)

        upwindc2D!(
            P,
            P_ex,
            vxc,
            vyc,
            NC,
            Δt,
            Δx,
            Δy,
        )

        # For negative vy and |CFLy| = 1:
        #
        #     Pᵢⱼⁿ⁺¹ = Pᵢ,ⱼ₊₁ⁿ
        #
        @test P[:,1:end-1] ≈ P0[:,2:end] atol=1e-12 rtol=1e-12

    end

end

@testset "Semi-Lagrangian linear-field translation" begin

    # ================================================================ #
    # Linear scalar field
    #
    #     P(x,y) = a + b*x + c*y
    #
    # advected by a constant velocity field.
    #
    # For constant velocity, the exact solution after one timestep is
    #
    #     Pⁿ⁺¹(x,y) =
    #         Pⁿ(x - vx*Δt, y - vy*Δt)
    #
    # Cubic interpolation should reproduce a linear function exactly.
    # ================================================================ #

    NC = (
        x = 8,
        y = 6,
    )

    Lx = 1.0
    Ly = 1.0

    Δx = Lx / NC.x
    Δy = Ly / NC.y

    # ------------------------------------------------------------ #
    # Coordinates
    # ------------------------------------------------------------ #

    xc = range(
        Δx / 2,
        Lx - Δx / 2;
        length=NC.x,
    )

    yc = range(
        Δy / 2,
        Ly - Δy / 2;
        length=NC.y,
    )

    # Extended centroid coordinates including one ghost layer
    xce = range(
        -Δx / 2,
        Lx + Δx / 2;
        length=NC.x + 2,
    )

    yce = range(
        -Δy / 2,
        Ly + Δy / 2;
        length=NC.y + 2,
    )

    x2d = repeat(collect(xc), 1, NC.y)
    y2d = repeat(collect(yc)', NC.x, 1)

    x = (
        c   = xc,
        ce  = xce,
        c2d = x2d,
    )

    y = (
        c   = yc,
        ce  = yce,
        c2d = y2d,
    )


    # ------------------------------------------------------------ #
    # Linear scalar field
    # ------------------------------------------------------------ #

    a = 2.0
    b = 3.0
    c = -1.5

    P = @. a + b * x2d + c * y2d

    # The extended field must represent the same linear function
    # including the ghost nodes.
    xe2d = repeat(collect(xce), 1, NC.y + 2)
    ye2d = repeat(collect(yce)', NC.x + 2, 1)

    P_ex = @. a + b * xe2d + c * ye2d


    # ------------------------------------------------------------ #
    # Constant velocity
    # ------------------------------------------------------------ #

    vx = 0.20
    vy = -0.10

    vxc = fill(vx, NC.x, NC.y)
    vyc = fill(vy, NC.x, NC.y)

    # Choose a sufficiently small timestep so all departure points
    # remain inside the extended scalar domain.
    Δt = 0.20

    # Exercise the branch in which no previous velocity field is
    # supplied. The routine then uses the current velocity directly.
    vxo = Float64[]
    vyo = Float64[]


    # ------------------------------------------------------------ #
    # Exact translated solution
    # ------------------------------------------------------------ #

    xd = @. x2d - vx * Δt
    yd = @. y2d - vy * Δt

    Pexact = @. a + b * xd + c * yd


    # ------------------------------------------------------------ #
    # Semi-Lagrangian update
    # ------------------------------------------------------------ #

    semilagc2D!(
        P,
        P_ex,
        vxc,
        vyc,
        vxo,
        vyo,
        x,
        y,
        Δt,
    )

    @test P ≈ Pexact atol=1e-12 rtol=1e-12

    # The interior of the extended field should be synchronized
    # with the newly advected scalar field.
    @test P_ex[2:end-1,2:end-1] ≈ Pexact atol=1e-12 rtol=1e-12

end
using GeoModBox.MomentumEquation.TwoD

@testset "2-D Stokes equation" begin

    # ================================================================ #
    # Pure shear, constant viscosity
    # ================================================================ #

    NC = (
        x = 6,
        y = 4,
    )

    NV = (
        x = NC.x + 1,
        y = NC.y + 1,
    )

    Lx = 2.0
    Ly = 1.0

    Δ = (
        x = Lx / NC.x,
        y = Ly / NC.y,
    )

    xmin = -Lx / 2
    xmax =  Lx / 2

    ymin = -Ly / 2
    ymax =  Ly / 2

    xv = range(xmin, xmax; length=NV.x)
    yv = range(ymin, ymax; length=NV.y)

    xc = range(
        xmin + Δ.x / 2,
        xmax - Δ.x / 2;
        length=NC.x,
    )

    yc = range(
        ymin + Δ.y / 2,
        ymax - Δ.y / 2;
        length=NC.y,
    )

    εbg = 2.0
    η   = 3.0
    g   = 0.0

    # ------------------------------------------------------------ #
    # Pure-shear boundary conditions
    # ------------------------------------------------------------ #

    BC = (
        type = (
            W = :ps,
            E = :ps,
            S = :ps,
            N = :ps,
        ),

        val = (
            # Tangential boundary values
            #
            # vx along S/N -> NV.x values
            # vy along W/E -> NV.y values
            W = zeros(NV.y),
            E = zeros(NV.y),
            S = zeros(NV.x),
            N = zeros(NV.x),

            # Normal velocities
            #
            # vx along W/E -> NC.y values
            # vy along S/N -> NC.x values
            vxW = fill(εbg * xmin, NC.y),
            vxE = fill(εbg * xmax, NC.y),

            vyS = fill(-εbg * ymin, NC.x),
            vyN = fill(-εbg * ymax, NC.x),
        ),
    )

    # ------------------------------------------------------------ #
    # Global numbering
    # ------------------------------------------------------------ #

    off = [
        NV.x * NC.y,
        NV.x * NC.y + NC.x * NV.y,
    ]

    Num = (
        Vx = reshape(
            1:NV.x * NC.y,
            NV.x,
            NC.y,
        ),

        Vy = reshape(
            off[1] + 1 : off[1] + NC.x * NV.y,
            NC.x,
            NV.y,
        ),

        Pt = reshape(
            off[2] + 1 : off[2] + NC.x * NC.y,
            NC.x,
            NC.y,
        ),
    )

    # ================================================================ #
    # Residual of the analytical pure-shear solution
    # ================================================================ #

    @testset "Pure shear residual" begin

        D = (
            # vx: x-vertices × extended y-centroids
            vx = zeros(NV.x, NC.y + 2),

            # vy: extended x-centroids × y-vertices
            vy = zeros(NC.x + 2, NV.y),

            Pt = zeros(NC.x, NC.y),
            ρ  = zeros(NC.x, NC.y),
        )

        # Analytical velocity field
        for i = 1:NV.x
            D.vx[i,2:end-1] .= εbg * xv[i]
        end

        for j = 1:NV.y
            D.vy[2:end-1,j] .= -εbg * yv[j]
        end

        divV = zeros(NC.x, NC.y)

        ε = (
            xx = zeros(NC.x, NC.y),
            yy = zeros(NC.x, NC.y),
            xy = zeros(NV.x, NV.y),
        )

        τ = (
            xx = zeros(NC.x, NC.y),
            yy = zeros(NC.x, NC.y),
            xy = zeros(NV.x, NV.y),
        )

        Fm = (
            x = zeros(NV.x, NC.y),
            y = zeros(NC.x, NV.y),
        )

        FPt = zeros(NC.x, NC.y)

        Residuals2Dc!(
            D,
            BC,
            ε,
            τ,
            divV,
            Δ,
            η,
            g,
            Fm,
            FPt,
        )

        # -------------------------------------------------------- #
        # Incompressibility
        # -------------------------------------------------------- #

        @test divV ≈ zeros(NC.x, NC.y) atol=1e-12 rtol=1e-12
        @test FPt  ≈ zeros(NC.x, NC.y) atol=1e-12 rtol=1e-12

        # -------------------------------------------------------- #
        # Strain rates
        # -------------------------------------------------------- #

        @test ε.xx ≈ fill( εbg, NC.x, NC.y) atol=1e-12 rtol=1e-12
        @test ε.yy ≈ fill(-εbg, NC.x, NC.y) atol=1e-12 rtol=1e-12
        @test ε.xy ≈ zeros(NV.x, NV.y) atol=1e-12 rtol=1e-12

        # -------------------------------------------------------- #
        # Stresses
        #
        # τij = 2 η εij
        # -------------------------------------------------------- #

        @test τ.xx ≈ fill( 2 * η * εbg, NC.x, NC.y) atol=1e-12 rtol=1e-12
        @test τ.yy ≈ fill(-2 * η * εbg, NC.x, NC.y) atol=1e-12 rtol=1e-12
        @test τ.xy ≈ zeros(NV.x, NV.y) atol=1e-12 rtol=1e-12

        # -------------------------------------------------------- #
        # Momentum residuals
        # -------------------------------------------------------- #

        @test Fm.x ≈ zeros(NV.x, NC.y) atol=1e-12 rtol=1e-12
        @test Fm.y ≈ zeros(NC.x, NV.y) atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Direct Stokes solve
    # ================================================================ #

    @testset "Pure shear direct solution" begin

        # ------------------------------------------------------------ #
        # Assemble system
        # ------------------------------------------------------------ #

        ρ = zeros(NC.x, NC.y)

        K = Assemblyc(
            NC,
            NV,
            Δ,
            η,
            BC,
            Num,
        )

        rhs = updaterhsc(
            NC,
            NV,
            Δ,
            η,
            ρ,
            g,
            BC,
            Num,
        )

        χ = K \ rhs

        # ------------------------------------------------------------ #
        # Extract numerical solution
        # ------------------------------------------------------------ #

        vx = χ[Num.Vx]
        vy = χ[Num.Vy]
        Pt = χ[Num.Pt]

        # ------------------------------------------------------------ #
        # Analytical solution at staggered nodes
        # ------------------------------------------------------------ #

        vx_exact = zeros(NV.x, NC.y)

        for i = 1:NV.x
            vx_exact[i,:] .= εbg * xv[i]
        end

        vy_exact = zeros(NC.x, NV.y)

        for j = 1:NV.y
            vy_exact[:,j] .= -εbg * yv[j]
        end

        # ------------------------------------------------------------ #
        # Velocity
        # ------------------------------------------------------------ #

        @test vx ≈ vx_exact atol=1e-11 rtol=1e-11
        @test vy ≈ vy_exact atol=1e-11 rtol=1e-11

        # One pressure DOF is fixed to remove the pressure nullspace.
        # Pure shear requires no pressure gradient, so all pressures
        # should consequently be approximately zero.
        @test Pt ≈ zeros(NC.x, NC.y) atol=1e-11 rtol=1e-11

    end


    # ================================================================ #
    # Boundary-condition validation
    # ================================================================ #

    @testset "Boundary-condition validation" begin

        BC_invalid = (
            type = (
                W = :freeslip,
                E = :freeslip,
                S = :unsupported,
                N = :freeslip,
            ),
        )

        @test_throws ArgumentError Assemblyc(
            NC,
            NV,
            Δ,
            η,
            BC_invalid,
            Num,
        )

    end

    @testset "Variable-viscosity solver: constant-viscosity limit" begin

        # ================================================================ #
        # Constant viscosity fields
        # ================================================================ #

        ηc = fill(η, NC.x, NC.y)
        ηv = fill(η, NV.x, NV.y)

        ρ = zeros(NC.x, NC.y)

        # ================================================================ #
        # Assemble and solve variable-viscosity system
        # ================================================================ #

        K = Assembly(
            NC,
            NV,
            Δ,
            ηc,
            ηv,
            BC,
            Num,
        )

        rhs = updaterhs(
            NC,
            NV,
            Δ,
            ηc,
            ηv,
            ρ,
            g,
            BC,
            Num,
        )

        χ = K \ rhs

        vx = χ[Num.Vx]
        vy = χ[Num.Vy]
        Pt = χ[Num.Pt]

        # ================================================================ #
        # Analytical pure-shear solution
        # ================================================================ #

        vx_exact = zeros(NV.x, NC.y)

        for i = 1:NV.x
            vx_exact[i,:] .= εbg * xv[i]
        end

        vy_exact = zeros(NC.x, NV.y)

        for j = 1:NV.y
            vy_exact[:,j] .= -εbg * yv[j]
        end

        @test vx ≈ vx_exact atol=1e-11 rtol=1e-11
        @test vy ≈ vy_exact atol=1e-11 rtol=1e-11
        @test Pt ≈ zeros(NC.x, NC.y) atol=1e-11 rtol=1e-11


        # ================================================================ #
        # Compare directly with constant-viscosity solver
        # ================================================================ #

        Kc = Assemblyc(
            NC,
            NV,
            Δ,
            η,
            BC,
            Num,
        )

        rhsc = updaterhsc(
            NC,
            NV,
            Δ,
            η,
            ρ,
            g,
            BC,
            Num,
        )

        χc = Kc \ rhsc

        vx_c = χc[Num.Vx]
        vy_c = χc[Num.Vy]
        Pt_c = χc[Num.Pt]

        @test vx ≈ vx_c atol=1e-11 rtol=1e-11
        @test vy ≈ vy_c atol=1e-11 rtol=1e-11
        @test Pt ≈ Pt_c atol=1e-11 rtol=1e-11


        # ================================================================ #
        # Variable-viscosity residual operator
        # ================================================================ #

        D = (
            vx = zeros(NV.x, NC.y + 2),
            vy = zeros(NC.x + 2, NV.y),
            Pt = copy(Pt),
            ρ  = copy(ρ),
        )

        D.vx[:,2:end-1] .= vx
        D.vy[2:end-1,:] .= vy

        divV = zeros(NC.x, NC.y)

        ε = (
            xx = zeros(NC.x, NC.y),
            yy = zeros(NC.x, NC.y),
            xy = zeros(NV.x, NV.y),
        )

        τ = (
            xx = zeros(NC.x, NC.y),
            yy = zeros(NC.x, NC.y),
            xy = zeros(NV.x, NV.y),
        )

        Fm = (
            x = zeros(NV.x, NC.y),
            y = zeros(NC.x, NV.y),
        )

        FPt = zeros(NC.x, NC.y)

        Residuals2D!(
            D,
            BC,
            ε,
            τ,
            divV,
            Δ,
            ηc,
            ηv,
            g,
            Fm,
            FPt,
        )

        @test divV ≈ zeros(NC.x, NC.y) atol=1e-11 rtol=1e-11
        @test FPt  ≈ zeros(NC.x, NC.y) atol=1e-11 rtol=1e-11
        @test Fm.x ≈ zeros(NV.x, NC.y) atol=1e-11 rtol=1e-11
        @test Fm.y ≈ zeros(NC.x, NV.y) atol=1e-11 rtol=1e-11

    end

    @testset "Hydrostatic equilibrium" begin

        # ================================================================ #
        # Constant density in a gravitational field
        #
        #     vx = vy = 0
        #
        # and vertical momentum balance requires
        #
        #     dP/dy = -ρg
        #
        # The pressure DOF Pt[1,1] is fixed to zero by the Stokes matrix,
        # so the analytical pressure is referenced to the first row of
        # cell centroids.
        # ================================================================ #

        ρ0 = 4.0
        gh = 3.0

        ρ = fill(ρ0, NC.x, NC.y)

        # ------------------------------------------------------------ #
        # Free-slip boundaries on all sides
        # ------------------------------------------------------------ #

        BC_hydro = (
            type = (
                W = :freeslip,
                E = :freeslip,
                S = :freeslip,
                N = :freeslip,
            ),

            val = (
                W = zeros(NV.y),
                E = zeros(NV.y),
                S = zeros(NV.x),
                N = zeros(NV.x),

                vxW = zeros(NC.y),
                vxE = zeros(NC.y),

                vyS = zeros(NC.x),
                vyN = zeros(NC.x),
            ),
        )


        # ================================================================ #
        # Constant-viscosity solver
        # ================================================================ #

        K = Assemblyc(
            NC,
            NV,
            Δ,
            η,
            BC_hydro,
            Num,
        )

        rhs = updaterhsc(
            NC,
            NV,
            Δ,
            η,
            ρ,
            gh,
            BC_hydro,
            Num,
        )

        χ = K \ rhs

        vx = χ[Num.Vx]
        vy = χ[Num.Vy]
        Pt = χ[Num.Pt]

        # No flow in hydrostatic equilibrium
        @test vx ≈ zeros(NV.x, NC.y) atol=1e-11 rtol=1e-11
        @test vy ≈ zeros(NC.x, NV.y) atol=1e-11 rtol=1e-11


        # ================================================================ #
        # Analytical hydrostatic pressure
        # ================================================================ #

        Pt_exact = zeros(NC.x, NC.y)

        for j = 1:NC.y
            Pt_exact[:,j] .= -ρ0 * gh * (yc[j] - yc[1])
        end

        @test Pt ≈ Pt_exact atol=1e-11 rtol=1e-11


        # ================================================================ #
        # Residual of analytical hydrostatic state
        # ================================================================ #

        D = (
            vx = zeros(NV.x, NC.y + 2),
            vy = zeros(NC.x + 2, NV.y),
            Pt = copy(Pt_exact),
            ρ  = copy(ρ),
        )

        divV = zeros(NC.x, NC.y)

        ε = (
            xx = zeros(NC.x, NC.y),
            yy = zeros(NC.x, NC.y),
            xy = zeros(NV.x, NV.y),
        )

        τ = (
            xx = zeros(NC.x, NC.y),
            yy = zeros(NC.x, NC.y),
            xy = zeros(NV.x, NV.y),
        )

        Fm = (
            x = zeros(NV.x, NC.y),
            y = zeros(NC.x, NV.y),
        )

        FPt = zeros(NC.x, NC.y)

        Residuals2Dc!(
            D,
            BC_hydro,
            ε,
            τ,
            divV,
            Δ,
            η,
            gh,
            Fm,
            FPt,
        )

        @test divV ≈ zeros(NC.x, NC.y) atol=1e-12 rtol=1e-12
        @test Fm.x ≈ zeros(NV.x, NC.y) atol=1e-12 rtol=1e-12
        @test Fm.y ≈ zeros(NC.x, NV.y) atol=1e-12 rtol=1e-12


        # ================================================================ #
        # Variable-viscosity solver
        #
        # Hydrostatic equilibrium should be independent of viscosity.
        # ================================================================ #

        ηc = fill(η, NC.x, NC.y)
        ηv = fill(η, NV.x, NV.y)

        K = Assembly(
            NC,
            NV,
            Δ,
            ηc,
            ηv,
            BC_hydro,
            Num,
        )

        rhs = updaterhs(
            NC,
            NV,
            Δ,
            ηc,
            ηv,
            ρ,
            gh,
            BC_hydro,
            Num,
        )

        χ = K \ rhs

        vx_var = χ[Num.Vx]
        vy_var = χ[Num.Vy]
        Pt_var = χ[Num.Pt]

        @test vx_var ≈ zeros(NV.x, NC.y) atol=1e-11 rtol=1e-11
        @test vy_var ≈ zeros(NC.x, NV.y) atol=1e-11 rtol=1e-11
        @test Pt_var ≈ Pt_exact atol=1e-11 rtol=1e-11

    end

end
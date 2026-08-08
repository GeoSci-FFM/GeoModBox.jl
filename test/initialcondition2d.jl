using GeoModBox.InitialCondition

@testset "2-D initial conditions" begin

    # ================================================================ #
    # Common grid
    # ================================================================ #

    NC = (
        x = 8,
        y = 6,
    )

    NV = (
        x = NC.x + 1,
        y = NC.y + 1,
    )

    M = (
        xmin = 0.0,
        xmax = 2.0,
        ymin = 0.0,
        ymax = 1.0,
    )

    Δ = (
        x = (M.xmax - M.xmin) / NC.x,
        y = (M.ymax - M.ymin) / NC.y,
    )

    # ------------------------------------------------------------ #
    # Coordinates
    # ------------------------------------------------------------ #

    xc = range(
        M.xmin + Δ.x / 2,
        M.xmax - Δ.x / 2;
        length=NC.x,
    )

    yc = range(
        M.ymin + Δ.y / 2,
        M.ymax - Δ.y / 2;
        length=NC.y,
    )

    xv = range(
        M.xmin,
        M.xmax;
        length=NV.x,
    )

    yv = range(
        M.ymin,
        M.ymax;
        length=NV.y,
    )

    xce = range(
        M.xmin - Δ.x / 2,
        M.xmax + Δ.x / 2;
        length=NC.x + 2,
    )

    yce = range(
        M.ymin - Δ.y / 2,
        M.ymax + Δ.y / 2;
        length=NC.y + 2,
    )

    # Extended staggered-grid coordinate arrays
    x_vx2d = repeat(collect(xv), 1, NC.y + 2)
    y_vx2d = repeat(collect(yce)', NV.x, 1)

    x_vy2d = repeat(collect(xce), 1, NV.y)
    y_vy2d = repeat(collect(yv)', NC.x + 2, 1)

    # Vertex-coordinate arrays
    x_v2d = repeat(collect(xv), 1, NV.y)
    y_v2d = repeat(collect(yv)', NV.x, 1)

    # Cell-center coordinate arrays
    x_c2d = repeat(collect(xc), 1, NC.y)
    y_c2d = repeat(collect(yc)', NC.x, 1)

    x = (
        c     = xc,
        v     = xv,
        ce    = xce,
        c2d   = x_c2d,
        v2d   = x_v2d,
        vx2d  = x_vx2d,
        vy2d  = x_vy2d,
    )

    y = (
        c     = yc,
        v     = yv,
        ce    = yce,
        c2d   = y_c2d,
        v2d   = y_v2d,
        vx2d  = y_vx2d,
        vy2d  = y_vy2d,
    )


    # ================================================================ #
    # Constant temperature
    # ================================================================ #

    @testset "Constant temperature" begin

        Tb = 750.0

        D = (
            T       = zeros(NC.x, NC.y),
            T_ex    = zeros(NC.x + 2, NC.y + 2),
            T0      = zeros(NC.x, NC.y),
            Told_ex = zeros(NC.x + 2, NC.y + 2),
            T_exD0  = zeros(NC.x + 2, NC.y + 2),
            T_ex0   = zeros(NC.x + 2, NC.y + 2),
        )

        IniTemperature!(
            :const,
            M,
            NC,
            D,
            x,
            y;
            Tb=Tb,
        )

        @test D.T_ex ≈ fill(Tb, NC.x + 2, NC.y + 2)
        @test D.T    ≈ fill(Tb, NC.x, NC.y)

        # Auxiliary fields should be synchronized as well
        @test D.T0      ≈ D.T
        @test D.Told_ex ≈ D.T_ex
        @test D.T_exD0  ≈ D.T_ex
        @test D.T_ex0   ≈ D.T_ex

    end


    # ================================================================ #
    # Linear temperature profile
    # ================================================================ #

    @testset "Linear temperature profile" begin

        Tbot = 1200.0
        Ttop = 300.0

        D = (
            T    = zeros(NC.x, NC.y),
            T_ex = zeros(NC.x + 2, NC.y + 2),
        )

        IniTemperature!(
            :linear,
            M,
            NC,
            D,
            x,
            y;
            Tb=Tbot,
            Ta=Ttop,
        )

        H = M.ymax - M.ymin

        Tgrad = (Tbot - Ttop) / H

        Texact = zeros(NC.x + 2, NC.y + 2)

        for i = 1:NC.x+2
            for j = 1:NC.y+2
                Texact[i,j] = -Tgrad * y.ce[j] + Ttop
            end
        end

        @test D.T_ex ≈ Texact atol=1e-12 rtol=1e-12
        @test D.T ≈ Texact[2:end-1,2:end-1] atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Blankenbach benchmark temperature
    # ================================================================ #

    @testset "Blankenbach temperature profile" begin

        Tbot = 1.0
        Ttop = 0.0

        D = (
            T    = zeros(NC.x, NC.y),
            T_ex = zeros(NC.x + 2, NC.y + 2),
        )

        IniTemperature!(
            :blankenbach,
            M,
            NC,
            D,
            x,
            y;
            Tb=Tbot,
            Ta=Ttop,
        )

        H = M.ymax - M.ymin
        L = M.xmax - M.xmin
        A = 0.01

        Texact = zeros(NC.x + 2, NC.y + 2)

        for i = 1:NC.x+2
            for j = 1:NC.y+2

                xn = (x.ce[i] - M.xmin) / L
                zn = (y.ce[j] - M.ymin) / H

                Tlinear =
                    Tbot + (Ttop - Tbot) * zn

                Texact[i,j] =
                    Tlinear +
                    A * (Tbot - Ttop) *
                    cos(π * xn) *
                    sin(π * zn)
            end
        end

        @test D.T_ex ≈ Texact atol=1e-12 rtol=1e-12
        @test D.T ≈ Texact[2:end-1,2:end-1] atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Pure-shear velocity initialization
    # ================================================================ #

    @testset "Pure-shear velocity boundaries" begin

        εbg = 2.0

        D = (
            vx = zeros(NV.x, NC.y + 2),
            vy = zeros(NC.x + 2, NV.y),
        )

        BC = (
            type = (
                W = :ps,
                E = :ps,
                S = :ps,
                N = :ps,
            ),

            val = (
                W   = zeros(NV.y),
                E   = zeros(NV.y),
                S   = zeros(NV.x),
                N   = zeros(NV.x),

                vxW = zeros(NC.y),
                vxE = zeros(NC.y),

                vyS = zeros(NC.x),
                vyN = zeros(NC.x),
            ),
        )

        D, BC = IniVelocity!(
            :PureShear,
            D,
            BC,
            NV,
            Δ,
            M,
            x,
            y;
            ε=εbg,
        )

        xc_domain = 0.5 * (M.xmin + M.xmax)
        yc_domain = 0.5 * (M.ymin + M.ymax)

        # ------------------------------------------------------------ #
        # Expected boundary values
        #
        # vx = -ε (x - xc)
        # vy =  ε (y - yc)
        # ------------------------------------------------------------ #

        S_exact = @. -(x.vx2d[:,1] - xc_domain) * εbg
        N_exact = @. -(x.vx2d[:,end] - xc_domain) * εbg

        vxW_exact = @. -(x.vx2d[1,2:end-1] - xc_domain) * εbg
        vxE_exact = @. -(x.vx2d[end,2:end-1] - xc_domain) * εbg

        vyS_exact = @. (y.vy2d[2:end-1,1] - yc_domain) * εbg
        vyN_exact = @. (y.vy2d[2:end-1,end] - yc_domain) * εbg

        W_exact = @. (y.v2d[1,:] - yc_domain) * εbg
        E_exact = @. (y.v2d[end,:] - yc_domain) * εbg

        @test BC.val.S   ≈ S_exact
        @test BC.val.N   ≈ N_exact

        @test BC.val.vxW ≈ vxW_exact
        @test BC.val.vxE ≈ vxE_exact

        @test BC.val.vyS ≈ vyS_exact
        @test BC.val.vyN ≈ vyN_exact

        @test BC.val.W   ≈ W_exact
        @test BC.val.E   ≈ E_exact

        # ------------------------------------------------------------ #
        # Boundary arrays should also have been written into D
        # ------------------------------------------------------------ #

        @test D.vx[:,1]         ≈ BC.val.S
        @test D.vx[:,end]       ≈ BC.val.N
        @test D.vx[1,2:end-1]   ≈ BC.val.vxW
        @test D.vx[end,2:end-1] ≈ BC.val.vxE

        @test D.vy[2:end-1,1]   ≈ BC.val.vyS
        @test D.vy[2:end-1,end] ≈ BC.val.vyN
        @test D.vy[1,:]         ≈ BC.val.W
        @test D.vy[end,:]       ≈ BC.val.E

    end


    # ================================================================ #
    # Block phase initialization
    # ================================================================ #

    @testset "Block phase initialization" begin

        phase = (0, 1)

        D = (
            p    = zeros(Int, NC.x, NC.y),
            p_ex = zeros(Int, NC.x + 2, NC.y + 2),
        )

        IniPhase!(
            :block,
            D,
            M,
            x,
            y,
            NC;
            phase=phase,
        )

        L = M.xmax - M.xmin
        H = M.ymax - M.ymin

        xL = M.xmin + 2/5 * L
        xR = M.xmin + 3/5 * L

        yU = M.ymax - 0.3 * H
        yO = M.ymax - 0.1 * H

        pexact = fill(phase[1], NC.x, NC.y)

        for i = 1:NC.x
            for j = 1:NC.y

                if x.c[i] >= xL &&
                   x.c[i] <= xR &&
                   y.c[j] >= yU &&
                   y.c[j] <= yO

                    pexact[i,j] = phase[2]
                end
            end
        end

        @test D.p == pexact

        # Both phases should actually be present
        @test phase[1] in D.p
        @test phase[2] in D.p

        # Interior synchronization
        @test D.p_ex[2:end-1,2:end-1] == D.p

        # Constant extrapolation into ghost cells
        @test D.p_ex[1,:]   == D.p_ex[2,:]
        @test D.p_ex[end,:] == D.p_ex[end-1,:]
        @test D.p_ex[:,1]   == D.p_ex[:,2]
        @test D.p_ex[:,end] == D.p_ex[:,end-1]

    end


    # ================================================================ #
    # Unsupported initial-condition types
    # ================================================================ #

    @testset "Unsupported initial-condition types" begin

        Dtemp = (
            T    = zeros(NC.x, NC.y),
            T_ex = zeros(NC.x + 2, NC.y + 2),
        )

        @test_throws ArgumentError IniTemperature!(
            :unsupported,
            M,
            NC,
            Dtemp,
            x,
            y,
        )


        Dvel = (
            vx = zeros(NV.x, NC.y + 2),
            vy = zeros(NC.x + 2, NV.y),
        )

        BC = (
            type = (
                W = :ps,
                E = :ps,
                S = :ps,
                N = :ps,
            ),

            val = (
                W   = zeros(NV.y),
                E   = zeros(NV.y),
                S   = zeros(NV.x),
                N   = zeros(NV.x),

                vxW = zeros(NC.y),
                vxE = zeros(NC.y),

                vyS = zeros(NC.x),
                vyN = zeros(NC.x),
            ),
        )

        @test_throws ArgumentError IniVelocity!(
            :unsupported,
            Dvel,
            BC,
            NV,
            Δ,
            M,
            x,
            y,
        )


        Dphase = (
            p    = zeros(Int, NC.x, NC.y),
            p_ex = zeros(Int, NC.x + 2, NC.y + 2),
        )

        @test_throws ArgumentError IniPhase!(
            :unsupported,
            Dphase,
            M,
            x,
            y,
            NC,
        )

    end

end
using GeoModBox.Tracers.TwoD

@testset "2-D tracers" begin

    # ================================================================ #
    # Grid setup
    # ================================================================ #

    NC = (
        x = 8,
        y = 6,
    )

    NV = (
        x = NC.x + 1,
        y = NC.y + 1,
    )

    Lx = 1.0
    Ly = 1.0

    Δ = (
        x = Lx / NC.x,
        y = Ly / NC.y,
    )

    # Cell centers
    xc = range(
        Δ.x / 2,
        Lx - Δ.x / 2;
        length=NC.x,
    )

    yc = range(
        Δ.y / 2,
        Ly - Δ.y / 2;
        length=NC.y,
    )

    # Vertices
    xv = range(
        0.0,
        Lx;
        length=NV.x,
    )

    yv = range(
        0.0,
        Ly;
        length=NV.y,
    )

    # Extended cell centers
    xce = range(
        -Δ.x / 2,
        Lx + Δ.x / 2;
        length=NC.x + 2,
    )

    yce = range(
        -Δ.y / 2,
        Ly + Δ.y / 2;
        length=NC.y + 2,
    )

    x = (
        c  = xc,
        v  = xv,
        ce = xce,
    )

    y = (
        c  = yc,
        v  = yv,
        ce = yce,
    )


    # ================================================================ #
    # Cell-centered field -> markers
    # ================================================================ #

    @testset "Cell-to-marker interpolation" begin

        # Linear field:
        #
        #     P(x,y) = a + b*x + c*y
        #
        # Bilinear interpolation should reproduce this exactly.

        a = 2.0
        b = 3.0
        c = -1.5

        xe2d = repeat(collect(xce), 1, NC.y + 2)
        ye2d = repeat(collect(yce)', NC.x + 2, 1)

        Prop = @. a + b * xe2d + c * ye2d

        Ma = Markers(
            [0.20, 0.45, 0.70],
            [0.25, 0.55, 0.75],
            [0, 0, 0],
        )

        for k = eachindex(Ma.x)

            Pm = FromCtoM(
                Prop,
                k,
                Ma,
                x,
                y,
                Δ,
                NC,
            )

            Pexact = a + b * Ma.x[k] + c * Ma.y[k]

            @test Pm ≈ Pexact atol=1e-12 rtol=1e-12
        end

    end


    # ================================================================ #
    # Staggered velocity -> markers
    # ================================================================ #

    @testset "Constant staggered velocity interpolation" begin

        vx0 = 0.25
        vy0 = -0.15

        # vx is defined on:
        #
        #     x.v × y.ce
        #
        Vx = fill(
            vx0,
            NC.x + 1,
            NC.y + 2,
        )

        # vy is defined on:
        #
        #     x.ce × y.v
        #
        Vy = fill(
            vy0,
            NC.x + 2,
            NC.y + 1,
        )

        Ma = Markers(
            [0.20, 0.45, 0.70],
            [0.25, 0.55, 0.75],
            [0, 0, 0],
        )

        for k = eachindex(Ma.x)

            vxm = VxFromVxNodes(
                Vx,
                k,
                Ma,
                x,
                y,
                Δ,
                NC,
                0,
            )

            vym = VyFromVyNodes(
                Vy,
                k,
                Ma,
                x,
                y,
                Δ,
                NC,
                0,
            )

            @test vxm ≈ vx0 atol=1e-12 rtol=1e-12
            @test vym ≈ vy0 atol=1e-12 rtol=1e-12
        end

    end


    # ================================================================ #
    # RK4 tracer advection
    # ================================================================ #

    @testset "RK4 constant-velocity translation" begin

        vx0 = 0.20
        vy0 = -0.10

        Vx = fill(
            vx0,
            NC.x + 1,
            NC.y + 2,
        )

        Vy = fill(
            vy0,
            NC.x + 2,
            NC.y + 1,
        )

        D = (
            vx = Vx,
            vy = Vy,
        )

        Ma = Markers(
            [0.25, 0.45, 0.65],
            [0.35, 0.55, 0.75],
            [0, 0, 0],
        )

        x0 = copy(Ma.x)
        y0 = copy(Ma.y)

        nmark = length(Ma.x)

        dt = 0.20

        rkw = (1.0 / 6.0) .* [1.0, 2.0, 2.0, 1.0]
        rkv = (1.0 / 2.0) .* [1.0, 1.0, 2.0, 2.0]

        AdvectTracer2D(
            Ma,
            nmark,
            D,
            x,
            y,
            dt,
            Δ,
            NC,
            rkw,
            rkv;
            style=1,
        )

        xexact = @. x0 + vx0 * dt
        yexact = @. y0 + vy0 * dt

        @test Ma.x ≈ xexact atol=1e-12 rtol=1e-12
        @test Ma.y ≈ yexact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Inactive markers
    # ================================================================ #

    @testset "Inactive markers remain fixed" begin

        vx0 = 0.20
        vy0 = -0.10

        D = (
            vx = fill(vx0, NC.x + 1, NC.y + 2),
            vy = fill(vy0, NC.x + 2, NC.y + 1),
        )

        # Middle marker is inactive
        Ma = Markers(
            [0.25, 0.45, 0.65],
            [0.35, 0.55, 0.75],
            [0, -1, 0],
        )

        x0 = copy(Ma.x)
        y0 = copy(Ma.y)

        dt = 0.20

        rkw = (1.0 / 6.0) .* [1.0, 2.0, 2.0, 1.0]
        rkv = (1.0 / 2.0) .* [1.0, 1.0, 2.0, 2.0]

        AdvectTracer2D(
            Ma,
            length(Ma.x),
            D,
            x,
            y,
            dt,
            Δ,
            NC,
            rkw,
            rkv;
            style=1,
        )

        @test Ma.x[2] == x0[2]
        @test Ma.y[2] == y0[2]

        @test Ma.x[1] ≈ x0[1] + vx0 * dt
        @test Ma.y[1] ≈ y0[1] + vy0 * dt

        @test Ma.x[3] ≈ x0[3] + vx0 * dt
        @test Ma.y[3] ≈ y0[3] + vy0 * dt

    end

    @testset "Marker-to-grid interpolation" begin

        # ================================================================ #
        # Grid
        # ================================================================ #

        NC = (
            x = 8,
            y = 6,
        )

        NV = (
            x = NC.x + 1,
            y = NC.y + 1,
        )

        Lx = 1.0
        Ly = 1.0

        Δ = (
            x = Lx / NC.x,
            y = Ly / NC.y,
        )

        xc = range(
            Δ.x / 2,
            Lx - Δ.x / 2;
            length=NC.x,
        )

        yc = range(
            Δ.y / 2,
            Ly - Δ.y / 2;
            length=NC.y,
        )

        xv = range(
            0.0,
            Lx;
            length=NV.x,
        )

        yv = range(
            0.0,
            Ly;
            length=NV.y,
        )

        xce = range(
            -Δ.x / 2,
            Lx + Δ.x / 2;
            length=NC.x + 2,
        )

        yce = range(
            -Δ.y / 2,
            Ly + Δ.y / 2;
            length=NC.y + 2,
        )

        x = (
            c  = xc,
            v  = xv,
            ce = xce,
        )

        y = (
            c  = yc,
            v  = yv,
            ce = yce,
        )


        # ================================================================ #
        # Marker positions
        #
        # Place markers at every cell centroid. Repeat the complete marker
        # distribution once per Julia thread. This guarantees that
        #
        #     nmark ÷ nthreads()
        #
        # is an integer and that Markers2Cells / Markers2Vertices create
        # exactly one chunk per thread.
        # ================================================================ #

        xm0 = vec(repeat(collect(xc), 1, NC.y))
        ym0 = vec(repeat(collect(yc)', NC.x, 1))

        nt = Base.Threads.maxthreadid()

        xm = repeat(xm0, nt)
        ym = repeat(ym0, nt)

        nmark = length(xm)


        # ================================================================ #
        # Constant marker temperature -> cells
        # ================================================================ #

        @testset "Constant temperature to cells" begin

            T0 = 750.0

            Ma = TMarkers(
                copy(xm),
                copy(ym),
                fill(T0, nmark),
                zeros(Int64, nmark),
            )

            # Extended centroid field
            PC = fill(T0, NC.x + 2, NC.y + 2)

            weight = zeros(NC.x + 2, NC.y + 2)

            PC_th = [
                zeros(NC.x + 2, NC.y + 2)
                for _ = 1:nt
            ]

            weight_th = [
                zeros(NC.x + 2, NC.y + 2)
                for _ = 1:nt
            ]

            Markers2Cells(
                Ma,
                nmark,
                PC_th,
                PC,
                weight_th,
                weight,
                x,
                y,
                Δ,
                :thermal,
                0,
            )

            @test PC ≈ fill(T0, size(PC)) atol=1e-12 rtol=1e-12

        end


        # ================================================================ #
        # Constant marker temperature -> vertices
        # ================================================================ #

        @testset "Constant temperature to vertices" begin

            T0 = 750.0

            Ma = TMarkers(
                copy(xm),
                copy(ym),
                fill(T0, nmark),
                zeros(Int64, nmark),
            )

            PG = fill(T0, NV.x, NV.y)

            weight = zeros(NV.x, NV.y)

            PG_th = [
                zeros(NV.x, NV.y)
                for _ = 1:nt
            ]

            weight_th = [
                zeros(NV.x, NV.y)
                for _ = 1:nt
            ]

            Markers2Vertices(
                Ma,
                nmark,
                PG_th,
                PG,
                weight_th,
                weight,
                x,
                y,
                Δ,
                :thermal,
                0,
            )

            @test PG ≈ fill(T0, size(PG)) atol=1e-12 rtol=1e-12

        end


        # ================================================================ #
        # Constant phase-dependent material property
        #
        # If every active marker carries the same phase and that phase has
        # property P0, arithmetic, harmonic, and geometric averaging must
        # all reproduce P0 exactly.
        # ================================================================ #

        @testset "Constant phase property to cells" begin

            P0 = 10.0

            Ma = Markers(
                copy(xm),
                copy(ym),
                zeros(Int64, nmark),
            )

            # Phase 0 -> param2[1]
            param2 = [P0]

            for avgm in (:arith, :harm, :geom)

                PC = fill(P0, NC.x + 2, NC.y + 2)

                weight = zeros(NC.x + 2, NC.y + 2)

                PC_th = [
                    zeros(NC.x + 2, NC.y + 2)
                    for _ = 1:nt
                ]

                weight_th = [
                    zeros(NC.x + 2, NC.y + 2)
                    for _ = 1:nt
                ]

                Markers2Cells(
                    Ma,
                    nmark,
                    PC_th,
                    PC,
                    weight_th,
                    weight,
                    x,
                    y,
                    Δ,
                    :phase,
                    param2;
                    avgm=avgm,
                )

                @test PC ≈ fill(P0, size(PC)) atol=1e-12 rtol=1e-12

            end

        end


        # ================================================================ #
        # Constant phase-dependent material property -> vertices
        # ================================================================ #

        @testset "Constant phase property to vertices" begin

            P0 = 10.0

            Ma = Markers(
                copy(xm),
                copy(ym),
                zeros(Int64, nmark),
            )

            param2 = [P0]

            for avgm in (:arith, :harm, :geom)

                PG = fill(P0, NV.x, NV.y)

                weight = zeros(NV.x, NV.y)

                PG_th = [
                    zeros(NV.x, NV.y)
                    for _ = 1:nt
                ]

                weight_th = [
                    zeros(NV.x, NV.y)
                    for _ = 1:nt
                ]

                Markers2Vertices(
                    Ma,
                    nmark,
                    PG_th,
                    PG,
                    weight_th,
                    weight,
                    x,
                    y,
                    Δ,
                    :phase,
                    param2;
                    avgm=avgm,
                )

                @test PG ≈ fill(P0, size(PG)) atol=1e-12 rtol=1e-12

            end

        end

    end

end


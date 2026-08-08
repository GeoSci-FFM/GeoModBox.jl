using GeoModBox.HeatEquation.TwoD

@testset "2-D heat diffusion" begin

    @testset "Constant field with zero-flux boundaries" begin

        # ------------------------------------------------------------ #
        # Grid
        # ------------------------------------------------------------ #

        NC = (
            x = 8,
            y = 6,
        )

        Δx = 1.0 / NC.x
        Δy = 1.0 / NC.y

        κ = 1.0

        # Stable timestep for Forward Euler:
        #
        # Δt <= 1 / [2κ(1/Δx² + 1/Δy²)]
        #
        Δt = 0.25 / (κ * (1 / Δx^2 + 1 / Δy^2))

        Tini = 1000.0

        # ------------------------------------------------------------ #
        # Boundary conditions
        # ------------------------------------------------------------ #

        BC = (
            type = (
                W = :Neumann,
                E = :Neumann,
                S = :Neumann,
                N = :Neumann,
            ),
            val = (
                W = zeros(NC.y),
                E = zeros(NC.y),
                S = zeros(NC.x),
                N = zeros(NC.x),
            ),
        )

        # Global numbering used by the matrix-based solvers
        Num = (
            T = reshape(1:NC.x * NC.y, NC.x, NC.y),
        )

        ndof = NC.x * NC.y

        Texact = fill(Tini, NC.x, NC.y)


        # ============================================================ #
        # Forward Euler
        # ============================================================ #

        D = (
            T    = copy(Texact),
            T_ex = fill(Tini, NC.x + 2, NC.y + 2),
        )

        ForwardEuler2Dc!(
            D,
            κ,
            Δx,
            Δy,
            Δt,
            NC,
            BC,
        )

        @test D.T ≈ Texact atol=1e-12 rtol=1e-12


        # ============================================================ #
        # Backward Euler
        # ============================================================ #

        D = (
            T    = copy(Texact),
            T_ex = fill(Tini, NC.x + 2, NC.y + 2),
        )

        rhs = zeros(ndof)
        K   = zeros(ndof, ndof)

        BackwardEuler2Dc!(
            D,
            κ,
            Δx,
            Δy,
            Δt,
            NC,
            BC,
            rhs,
            K,
            Num,
        )

        @test D.T ≈ Texact atol=1e-12 rtol=1e-12


        # ============================================================ #
        # Crank-Nicolson
        # ============================================================ #

        D = (
            T    = copy(Texact),
            T_ex = fill(Tini, NC.x + 2, NC.y + 2),
        )

        rhs = zeros(ndof)
        K1  = zeros(ndof, ndof)
        K2  = zeros(ndof, ndof)

        CNA2Dc!(
            D,
            κ,
            Δx,
            Δy,
            Δt,
            NC,
            BC,
            rhs,
            K1,
            K2,
            Num,
        )

        @test D.T ≈ Texact atol=1e-12 rtol=1e-12


        # ============================================================ #
        # ADI
        # ============================================================ #

        D = (
            T    = copy(Texact),
            T_ex = fill(Tini, NC.x + 2, NC.y + 2),
        )

        ADI2Dc!(
            D,
            κ,
            Δx,
            Δy,
            Δt,
            NC,
            BC,
        )

        @test D.T ≈ Texact atol=1e-12 rtol=1e-12

    end

end

@testset "Uniform volumetric heating" begin

    # ================================================================ #
    # Uniform temperature + zero-flux boundaries + uniform Q
    #
    # Since ∇²T = 0,
    #
    #     Tⁿ⁺¹ = Tⁿ + Q Δt / (ρ cp)
    #
    # ================================================================ #

    NC = (
        x = 8,
        y = 6,
    )

    Δx = 1.0 / NC.x
    Δy = 1.0 / NC.y

    κ = 1.0

    Δt = 0.25 / (κ * (1 / Δx^2 + 1 / Δy^2))

    Tini = 1000.0

    ρ  = 3300.0
    cp = 1250.0

    Q0 = 2.5e6
    Q  = fill(Q0, NC.x, NC.y)

    BC = (
        type = (
            W = :Neumann,
            E = :Neumann,
            S = :Neumann,
            N = :Neumann,
        ),
        val = (
            W = zeros(NC.y),
            E = zeros(NC.y),
            S = zeros(NC.x),
            N = zeros(NC.x),
        ),
    )

    Num = (
        T = reshape(1:NC.x * NC.y, NC.x, NC.y),
    )

    ndof = NC.x * NC.y

    Texact = fill(
        Tini + Q0 * Δt / (ρ * cp),
        NC.x,
        NC.y,
    )


    # ============================================================ #
    # Forward Euler
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    ForwardEuler2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # Backward Euler
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    rhs = zeros(ndof)
    K   = zeros(ndof, ndof)

    BackwardEuler2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
        rhs,
        K,
        Num;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # Crank-Nicolson
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    rhs = zeros(ndof)
    K1  = zeros(ndof, ndof)
    K2  = zeros(ndof, ndof)

    CNA2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
        rhs,
        K1,
        K2,
        Num;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # ADI
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    ADI2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12

end

@testset "Uniform shear heating" begin

    # ================================================================ #
    # Uniform temperature + zero-flux boundaries + uniform Qₛ
    #
    #     Tⁿ⁺¹ = Tⁿ + Qₛ Δt / (ρ cp)
    #
    # ADI2Dc! is not included because it currently has no Qₛ keyword.
    # ================================================================ #

    NC = (
        x = 8,
        y = 6,
    )

    Δx = 1.0 / NC.x
    Δy = 1.0 / NC.y

    κ = 1.0

    Δt = 0.25 / (κ * (1 / Δx^2 + 1 / Δy^2))

    Tini = 1000.0

    ρ  = 3300.0
    cp = 1250.0

    Qs0 = 1.5e6
    Qₛ  = fill(Qs0, NC.x, NC.y)

    BC = (
        type = (
            W = :Neumann,
            E = :Neumann,
            S = :Neumann,
            N = :Neumann,
        ),
        val = (
            W = zeros(NC.y),
            E = zeros(NC.y),
            S = zeros(NC.x),
            N = zeros(NC.x),
        ),
    )

    Num = (
        T = reshape(1:NC.x * NC.y, NC.x, NC.y),
    )

    ndof = NC.x * NC.y

    Texact = fill(
        Tini + Qs0 * Δt / (ρ * cp),
        NC.x,
        NC.y,
    )


    # ============================================================ #
    # Forward Euler
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    ForwardEuler2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC;
        Qₛ=Qₛ,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # Backward Euler
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    rhs = zeros(ndof)
    K   = zeros(ndof, ndof)

    BackwardEuler2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
        rhs,
        K,
        Num;
        Qₛ=Qₛ,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # Crank-Nicolson
    # ============================================================ #

    D = (
        T    = fill(Tini, NC.x, NC.y),
        T_ex = fill(Tini, NC.x + 2, NC.y + 2),
    )

    rhs = zeros(ndof)
    K1  = zeros(ndof, ndof)
    K2  = zeros(ndof, ndof)

    CNA2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
        rhs,
        K1,
        K2,
        Num;
        Qₛ=Qₛ,
        ρ=ρ,
        cp=cp,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12

end

@testset "Linear steady-state profile with mixed boundaries" begin

    # ================================================================ #
    # Steady linear temperature profile
    #
    #     T(y) = TS + (TN - TS) * y / H
    #
    # with
    #
    #     South/North : Dirichlet
    #     West/East   : Neumann
    #
    # Since ∂²T/∂x² = ∂²T/∂y² = 0, the solution should remain unchanged.
    # ================================================================ #

    NC = (
        x = 8,
        y = 6,
    )

    Lx = 2.0
    Ly = 1.0

    Δx = Lx / NC.x
    Δy = Ly / NC.y

    κ = 1.0

    Δt = 0.25 / (κ * (1 / Δx^2 + 1 / Δy^2))

    TS = 300.0
    TN = 1300.0

    yc = range(Δy / 2, Ly - Δy / 2; length=NC.y)

    Texact = zeros(NC.x, NC.y)

    for j = 1:NC.y
        @views Texact[:,j] .= TS + (TN - TS) * yc[j] / Ly
    end

    BC = (
        type = (
            W = :Neumann,
            E = :Neumann,
            S = :Dirichlet,
            N = :Dirichlet,
        ),
        val = (
            W = zeros(NC.y),
            E = zeros(NC.y),
            S = fill(TS, NC.x),
            N = fill(TN, NC.x),
        ),
    )

    Num = (
        T = reshape(1:NC.x * NC.y, NC.x, NC.y),
    )

    ndof = NC.x * NC.y


    # ============================================================ #
    # Forward Euler
    # ============================================================ #

    D = (
        T    = copy(Texact),
        T_ex = zeros(NC.x + 2, NC.y + 2),
    )

    D.T_ex[2:end-1,2:end-1] .= D.T

    ForwardEuler2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # Backward Euler
    # ============================================================ #

    D = (
        T    = copy(Texact),
        T_ex = zeros(NC.x + 2, NC.y + 2),
    )

    rhs = zeros(ndof)
    K   = zeros(ndof, ndof)

    BackwardEuler2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
        rhs,
        K,
        Num,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # Crank-Nicolson
    # ============================================================ #

    D = (
        T    = copy(Texact),
        T_ex = zeros(NC.x + 2, NC.y + 2),
    )

    rhs = zeros(ndof)
    K1  = zeros(ndof, ndof)
    K2  = zeros(ndof, ndof)

    CNA2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
        rhs,
        K1,
        K2,
        Num,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12


    # ============================================================ #
    # ADI
    # ============================================================ #

    D = (
        T    = copy(Texact),
        T_ex = zeros(NC.x + 2, NC.y + 2),
    )

    ADI2Dc!(
        D,
        κ,
        Δx,
        Δy,
        Δt,
        NC,
        BC,
    )

    @test D.T ≈ Texact atol=1e-12 rtol=1e-12

end
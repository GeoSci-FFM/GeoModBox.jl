using GeoModBox.HeatEquation.OneD

@testset "1-D heat diffusion" begin

    # ------------------------------------------------------------------ #
    # Constant temperature field with zero-flux boundaries
    #
    # For
    #
    #     ∂T/∂x = 0
    #     Q     = 0
    #
    # a spatially constant temperature field must remain constant.
    # ------------------------------------------------------------------ #

    nc  = 16
    Δx  = 1.0 / nc
    κ   = 1.0

    # Stable timestep for Forward Euler
    Δt  = 0.25 * Δx^2 / κ

    Tini = 1000.0

    BC = (
        type = (
            W = :Neumann,
            E = :Neumann,
        ),
        val = (
            W = 0.0,
            E = 0.0,
        ),
    )

    # ================================================================ #
    # Forward Euler
    # ================================================================ #

    explicit = (
        T    = fill(Tini, nc),
        T_ex = fill(Tini, nc + 2),
    )

    ForwardEuler1Dc!(
        explicit,
        κ,
        Δx,
        Δt,
        nc,
        BC,
    )

    @test all(isapprox.(explicit.T, Tini; atol=1e-12, rtol=0.0))
    @test all(isapprox.(explicit.T_ex[2:end-1], Tini; atol=1e-12, rtol=0.0))

    # ================================================================ #
    # Backward Euler
    # ================================================================ #

    implicit = (
        T = fill(Tini, nc),
    )

    K   = zeros(nc, nc)
    rhs = zeros(nc)

    BackwardEuler1Dc!(
        implicit,
        κ,
        Δx,
        Δt,
        nc,
        BC,
        K,
        rhs,
    )

    @test all(isapprox.(implicit.T, Tini; atol=1e-12, rtol=0.0))

    # ================================================================ #
    # Crank-Nicolson
    # ================================================================ #

    cna = (
        T = fill(Tini, nc),
    )

    K1 = zeros(nc, nc)
    K2 = zeros(nc, nc)

    CNA1Dc!(
        cna,
        κ,
        Δx,
        Δt,
        nc,
        BC,
        K1,
        K2,
    )

    @test all(isapprox.(cna.T, Tini; atol=1e-12, rtol=0.0))

end

@testset "Uniform internal heating" begin

    nc  = 16
    Δx  = 1.0 / nc
    κ   = 1.0
    Δt  = 0.25 * Δx^2 / κ

    Tini = 1000.0

    ρ  = 3300.0
    cp = 1250.0
    Q  = fill(2.5e6, nc)

    BC = (
        type = (
            W = :Neumann,
            E = :Neumann,
        ),
        val = (
            W = 0.0,
            E = 0.0,
        ),
    )

    # Exact temperature after one timestep
    Texact = Tini + Q[1] * Δt / (ρ * cp)

    # ================================================================ #
    # Forward Euler
    # ================================================================ #

    explicit = (
        T    = fill(Tini, nc),
        T_ex = fill(Tini, nc + 2),
    )

    ForwardEuler1Dc!(
        explicit,
        κ,
        Δx,
        Δt,
        nc,
        BC;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test all(isapprox.(explicit.T, Texact; atol=1e-12, rtol=1e-12))

    # ================================================================ #
    # Backward Euler
    # ================================================================ #

    implicit = (
        T = fill(Tini, nc),
    )

    K   = zeros(nc, nc)
    rhs = zeros(nc)

    BackwardEuler1Dc!(
        implicit,
        κ,
        Δx,
        Δt,
        nc,
        BC,
        K,
        rhs;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test all(isapprox.(implicit.T, Texact; atol=1e-12, rtol=1e-12))

    # ================================================================ #
    # Crank-Nicolson
    # ================================================================ #

    cna = (
        T = fill(Tini, nc),
    )

    K1 = zeros(nc, nc)
    K2 = zeros(nc, nc)

    CNA1Dc!(
        cna,
        κ,
        Δx,
        Δt,
        nc,
        BC,
        K1,
        K2;
        Q=Q,
        ρ=ρ,
        cp=cp,
    )

    @test all(isapprox.(cna.T, Texact; atol=1e-12, rtol=1e-12))

end

@testset "Linear steady-state profile with Dirichlet boundaries" begin

    # ------------------------------------------------------------------ #
    # A linear temperature profile between two fixed boundary
    # temperatures is a steady-state solution of
    #
    #     ∂T/∂t = κ ∂²T/∂x²
    #
    # because ∂²T/∂x² = 0.
    #
    # The numerical solution should therefore remain unchanged.
    # ------------------------------------------------------------------ #

    nc = 16
    L  = 1.0
    Δx = L / nc
    κ  = 1.0

    # Stable timestep for Forward Euler
    Δt = 0.25 * Δx^2 / κ

    TW = 300.0
    TE = 1300.0

    # Cell-centroid coordinates
    x = range(Δx / 2, L - Δx / 2; length=nc)

    # Exact steady-state solution
    Texact = @. TW + (TE - TW) * x / L

    BC = (
        type = (
            W = :Dirichlet,
            E = :Dirichlet,
        ),
        val = (
            W = TW,
            E = TE,
        ),
    )

    # ================================================================ #
    # Forward Euler
    # ================================================================ #

    explicit = (
        T    = copy(Texact),
        T_ex = zeros(nc + 2),
    )

    explicit.T_ex[2:end-1] .= explicit.T

    ForwardEuler1Dc!(
        explicit,
        κ,
        Δx,
        Δt,
        nc,
        BC,
    )

    @test explicit.T ≈ Texact atol=1e-12 rtol=1e-12

    # ================================================================ #
    # Backward Euler
    # ================================================================ #

    implicit = (
        T = copy(Texact),
    )

    K   = zeros(nc, nc)
    rhs = zeros(nc)

    BackwardEuler1Dc!(
        implicit,
        κ,
        Δx,
        Δt,
        nc,
        BC,
        K,
        rhs,
    )

    @test implicit.T ≈ Texact atol=1e-12 rtol=1e-12

    # ================================================================ #
    # Crank-Nicolson
    # ================================================================ #

    cna = (
        T = copy(Texact),
    )

    K1 = zeros(nc, nc)
    K2 = zeros(nc, nc)

    CNA1Dc!(
        cna,
        κ,
        Δx,
        Δt,
        nc,
        BC,
        K1,
        K2,
    )

    @test cna.T ≈ Texact atol=1e-12 rtol=1e-12

end
using GeoModBox.Scaling

# ------------------------------------------------------------------ #
# Minimal mutable structures used only for testing ScaleParameters!
# ------------------------------------------------------------------ #

mutable struct TestScalingGeometry
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
end

mutable struct TestScalingSpacing
    x::Float64
    y::Float64
end

mutable struct TestScalingTime
    tmax::Float64
    Δc::Float64
    Δd::Float64
    Δ::Float64
end

mutable struct TestScalingPhysics
    κ::Float64
    η₀::Float64
    ΔT::Float64
    ρ₀::Float64
    cp::Float64
    Ttop::Float64
    Tbot::Float64
end

mutable struct TestScalingData
    Q::Vector{Float64}
end


@testset "Scaling" begin

    # ================================================================ #
    # Characteristic scaling constants
    # ================================================================ #

    @testset "Scaling constants" begin

        M = TestScalingGeometry(
            0.0,
            6.0,
            -2.0,
            0.0,
        )

        P = TestScalingPhysics(
            0.5,        # κ
            8.0,        # η₀
            1000.0,     # ΔT
            4.0,        # ρ₀
            5.0,        # cp
            273.15,     # Ttop
            1273.15,    # Tbot
        )

        S = ScalingConstants!(M, P)

        # Model height
        H = M.ymax - M.ymin

        @test S.hsc ≈ H

        # Diffusion time scale: H² / κ
        @test S.tsc ≈ H^2 / P.κ

        # Velocity scale: κ / H
        @test S.vsc ≈ P.κ / H

        # Stress scale: η₀ κ / H²
        @test S.τsc ≈ P.η₀ * P.κ / H^2

        # Temperature scale
        @test S.Tsc ≈ P.ΔT

        # Volumetric heat-production scale
        @test S.Qsc ≈ P.ΔT * P.κ * P.ρ₀ * P.cp / H^2

    end


    # ================================================================ #
    # Scaling of dimensional parameters
    # ================================================================ #

    @testset "Scale parameters" begin

        M = TestScalingGeometry(
            0.0,
            6.0,
            -2.0,
            0.0,
        )

        Δ = TestScalingSpacing(
            0.5,
            0.25,
        )

        T = TestScalingTime(
            16.0,       # tmax
            4.0,        # Δc
            2.0,        # Δd
            1.0,        # Δ
        )

        P = TestScalingPhysics(
            0.5,        # κ
            8.0,        # η₀
            1000.0,     # ΔT
            4.0,        # ρ₀
            5.0,        # cp
            273.15,     # Ttop
            1273.15,    # Tbot
        )

        S = ScalingConstants!(M, P)

        # Choose dimensional heat production equal to one and two
        # characteristic heat-production units.
        D = TestScalingData(
            [S.Qsc, 2.0 * S.Qsc],
        )

        # Save expected values before scaling in place
        xmin_expected = M.xmin / S.hsc
        xmax_expected = M.xmax / S.hsc
        ymin_expected = M.ymin / S.hsc
        ymax_expected = M.ymax / S.hsc

        Δx_expected = Δ.x / S.hsc
        Δy_expected = Δ.y / S.hsc

        tmax_expected = T.tmax / S.tsc
        Δc_expected   = T.Δc   / S.tsc
        Δd_expected   = T.Δd   / S.tsc
        Δt_expected   = T.Δ    / S.tsc

        Ttop_expected = (P.Ttop - 273.15) / S.Tsc
        Tbot_expected = (P.Tbot - 273.15) / S.Tsc

        ScaleParameters!(S, M, Δ, T, P, D)

        # ------------------------------------------------------------ #
        # Geometry
        # ------------------------------------------------------------ #

        @test M.xmin ≈ xmin_expected
        @test M.xmax ≈ xmax_expected
        @test M.ymin ≈ ymin_expected
        @test M.ymax ≈ ymax_expected

        @test Δ.x ≈ Δx_expected
        @test Δ.y ≈ Δy_expected

        # ------------------------------------------------------------ #
        # Time
        # ------------------------------------------------------------ #

        @test T.tmax ≈ tmax_expected
        @test T.Δc   ≈ Δc_expected
        @test T.Δd   ≈ Δd_expected
        @test T.Δ    ≈ Δt_expected

        # ------------------------------------------------------------ #
        # Temperature
        # ------------------------------------------------------------ #

        @test P.Ttop ≈ Ttop_expected
        @test P.Tbot ≈ Tbot_expected

        # For the chosen temperatures:
        #
        # Ttop = 273.15 K  -> 0
        # Tbot = 1273.15 K -> 1
        #
        @test P.Ttop ≈ 0.0
        @test P.Tbot ≈ 1.0

        # ------------------------------------------------------------ #
        # Heat production
        # ------------------------------------------------------------ #

        @test D.Q ≈ [1.0, 2.0]

    end

end
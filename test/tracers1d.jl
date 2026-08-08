using GeoModBox.Tracers.OneD

@testset "1-D tracers" begin

    # ================================================================ #
    # Centroids -> markers
    # ================================================================ #

    @testset "Center-to-marker interpolation" begin

        Δx   = 1.0
        xmin = 0.0

        xc = [0.5, 1.5, 2.5, 3.5]

        # Linear centroid field:
        #
        #     T(x) = 2 + 3x
        #
        Tc = @. 2.0 + 3.0 * xc

        # Include markers both inside and outside the centroid interval
        # to test interpolation and linear extrapolation.
        xm = [-0.5, 1.0, 2.0, 3.0, 4.5]

        Tm = zeros(length(xm))

        Itp1D_Centers2Markers!(
            Tm,
            xm,
            Tc,
            xc,
            Δx,
            xmin,
        )

        Texact = @. 2.0 + 3.0 * xm

        @test Tm ≈ Texact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Markers -> centroids
    # ================================================================ #

    @testset "Marker-to-center interpolation" begin

        Δx   = 1.0
        xmin = 0.0

        xc = [0.5, 1.5, 2.5, 3.5]

        # Put one marker exactly at every centroid.
        xm = copy(xc)

        # Linear marker field
        Tm = @. 2.0 + 3.0 * xm

        Tc = zeros(length(xc))

        Itp1D_Markers2Centers!(
            Tc,
            xc,
            Tm,
            xm,
            Δx,
            xmin,
        )

        Texact = @. 2.0 + 3.0 * xc

        @test Tc ≈ Texact atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Empty centroids retain their previous values
    # ================================================================ #

    @testset "Empty centers retain old values" begin

        Δx   = 1.0
        xmin = 0.0

        xc = [0.5, 1.5, 2.5, 3.5]

        # Only the first two centroids receive marker contributions.
        xm = [0.5, 1.5]
        Tm = [100.0, 200.0]

        # Existing values in cells without markers should survive.
        Tc = [10.0, 20.0, 30.0, 40.0]

        Itp1D_Markers2Centers!(
            Tc,
            xc,
            Tm,
            xm,
            Δx,
            xmin,
        )

        @test Tc[1] ≈ 100.0
        @test Tc[2] ≈ 200.0

        # No marker contributed to these cells.
        @test Tc[3] ≈ 30.0
        @test Tc[4] ≈ 40.0

    end


    # ================================================================ #
    # Constant marker field
    # ================================================================ #

    @testset "Constant marker field" begin

        Δx   = 1.0
        xmin = 0.0

        xc = [0.5, 1.5, 2.5, 3.5]

        # Markers located between the centroids.
        xm = [
            0.5,
            1.0,
            1.5,
            2.0,
            2.5,
            3.0,
            3.5,
        ]

        T0 = 750.0

        Tm = fill(T0, length(xm))
        Tc = zeros(length(xc))

        Itp1D_Markers2Centers!(
            Tc,
            xc,
            Tm,
            xm,
            Δx,
            xmin,
        )

        # A weighted average of a constant field must remain constant.
        @test Tc ≈ fill(T0, length(xc)) atol=1e-12 rtol=1e-12

    end


    # ================================================================ #
    # Input validation
    # ================================================================ #

    @testset "Input validation" begin

        xc = [0.5, 1.5, 2.5, 3.5]
        Tc = [10.0, 20.0, 30.0, 40.0]

        xm = [1.0, 2.0]
        Tm = zeros(2)

        xmin = 0.0

        # Tm and xm must have identical lengths
        @test_throws DimensionMismatch Itp1D_Centers2Markers!(
            zeros(3),
            xm,
            Tc,
            xc,
            1.0,
            xmin,
        )

        # Tc and xc must have identical lengths
        @test_throws DimensionMismatch Itp1D_Centers2Markers!(
            Tm,
            xm,
            Tc[1:3],
            xc,
            1.0,
            xmin,
        )

        # At least two centroids are required
        @test_throws ArgumentError Itp1D_Centers2Markers!(
            Tm,
            xm,
            [10.0],
            [0.5],
            1.0,
            xmin,
        )

        # Grid spacing must be positive
        @test_throws ArgumentError Itp1D_Centers2Markers!(
            Tm,
            xm,
            Tc,
            xc,
            0.0,
            xmin,
        )

        # Test the corresponding Δx protection of the reverse
        # interpolation routine as well.
        @test_throws ArgumentError Itp1D_Markers2Centers!(
            copy(Tc),
            xc,
            [10.0, 20.0],
            xm,
            -1.0,
            xmin,
        )

    end

end
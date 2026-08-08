using Test
using GeoModBox

function run_testfile(file, description)
    @info "Testing $description..."

    t = @elapsed include(file)

    @info "Finished $description" time=round(t; digits=2)
end

@testset "GeoModBox.jl" begin

    @info "Testing package loading..."

    @testset "Package loading" begin
        @test isdefined(GeoModBox, :HeatEquation)
        @test isdefined(GeoModBox, :AdvectionEquation)
        @test isdefined(GeoModBox, :MomentumEquation)
        @test isdefined(GeoModBox, :InitialCondition)
        @test isdefined(GeoModBox, :Tracers)
        @test isdefined(GeoModBox, :Scaling)
    end

    @info "Testing core structures..."

    @testset "Structures" begin
        M = Geometry()

        @test M isa Geometry
        @test M.xmax > M.xmin
        @test M.ymax > M.ymin

        T = TimeParameter()

        @test T.year ≈ 365.25 * 24 * 3600
        @test T.Δ == 0.0
        @test T.itmax > 0
    end

    run_testfile("heat1d.jl", "1-D heat equation")
    run_testfile("advection1d.jl", "1-D advection")
    run_testfile("tracers1d.jl", "1-D tracers")
    run_testfile("scaling.jl", "scaling")
    run_testfile("momentum1d.jl", "1-D momentum equation")

    run_testfile("initialcondition2d.jl", "2-D initial conditions")
    run_testfile("heat2d.jl", "2-D heat equation")
    run_testfile("advection2d.jl", "2-D advection")
    run_testfile("tracers2d.jl", "2-D tracers")
    run_testfile("momentum2d.jl", "2-D momentum equation")

end

# using Test
# using GeoModBox

# @testset "GeoModBox.jl" begin

#     @testset "Package loading" begin
#         @test isdefined(GeoModBox, :HeatEquation)
#         @test isdefined(GeoModBox, :AdvectionEquation)
#         @test isdefined(GeoModBox, :MomentumEquation)
#         @test isdefined(GeoModBox, :InitialCondition)
#         @test isdefined(GeoModBox, :Tracers)
#         @test isdefined(GeoModBox, :Scaling)
#     end

#     @testset "Structures" begin
#         M = Geometry()

#         @test M isa Geometry
#         @test M.xmax > M.xmin
#         @test M.ymax > M.ymin

#         T = TimeParameter()

#         @test T.year ≈ 365.25 * 24 * 3600
#         @test T.Δ == 0.0
#         @test T.itmax > 0
#     end

#     include("heat1d.jl")
#     include("advection1d.jl")
#     include("tracers1d.jl")
#     include("scaling.jl")
#     include("momentum1d.jl")
#     include("initialcondition2d.jl")
#     include("heat2d.jl")
#     include("advection2d.jl")
#     include("tracers2d.jl")
#     include("momentum2d.jl")

# end
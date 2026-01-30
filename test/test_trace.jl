using Test
using GeoCotrans
using FieldTracer
using Dates
using LinearAlgebra

@testset "trace_field_line" begin
    t = DateTime(2020, 1, 1)

    @testset "Basic tracing from dayside equator" begin
        # Start at L=4 on dayside equator (GSM coordinates)
        x, y, z = 4.0, 0.0, 0.0
        result = trace_field_line(x, y, z, t; rlim=10.0, r0=1.0)

        # Should trace to inner boundary
        @test result.status == :inner_boundary

        # Should have multiple points
        @test length(result) > 10

        # First point should be the starting position
        @test result.points[1] ≈ [x, y, z]

        # Last point should be at r ≈ r0
        @test result.r[end] ≈ 1.0 rtol=0.1

        # All intermediate radii should be > r0
        @test all(r -> r >= 0.99, result.r)
    end

    @testset "Tracing in both directions" begin
        x, y, z = 3.0, 0.0, 0.5

        # Trace parallel to B
        result_fwd = trace_field_line(x, y, z, t; dir=1, rlim=10.0, r0=1.0)

        # Trace anti-parallel to B
        result_bwd = trace_field_line(x, y, z, t; dir=-1, rlim=10.0, r0=1.0)

        # Both should reach inner boundary for dipole-like field
        @test result_fwd.status == :inner_boundary
        @test result_bwd.status == :inner_boundary

        # The endpoints should be in different hemispheres (approximately)
        @test sign(result_fwd.points[end][3]) != sign(result_bwd.points[end][3]) ||
              abs(result_fwd.points[end][3]) < 0.5 || abs(result_bwd.points[end][3]) < 0.5
    end

    @testset "Outer boundary reached" begin
        # Start far from Earth with small rlim
        x, y, z = -5.0, 0.0, 3.0
        result = trace_field_line(x, y, z, t; rlim=6.0, r0=1.0, dir=1)

        # Should eventually hit outer boundary or inner boundary
        @test result.status in (:outer_boundary, :inner_boundary, :lateral_boundary)
    end

    @testset "Vector input interface" begin
        pos = [3.0, 0.0, 0.0]
        result = trace_field_line(pos, t; rlim=10.0, r0=1.0)

        @test result.status == :inner_boundary
        @test length(result) > 10
    end

    @testset "CoordinateVector input interface" begin
        pos = GSM(3.0, 0.0, 0.0)
        result = trace_field_line(pos, t; rlim=10.0, r0=1.0)

        @test result.status == :inner_boundary
        @test length(result) > 10
    end

    @testset "Different coordinate systems" begin
        # Test GEO input
        pos_geo = GEO(3.0, 0.0, 0.0)
        result_geo = trace_field_line(pos_geo, t; rlim=10.0, r0=1.0)
        @test result_geo.status == :inner_boundary

        # Results should differ since the starting position is different in GSM
        pos_gsm = GSM(3.0, 0.0, 0.0)
        result_gsm = trace_field_line(pos_gsm, t; rlim=10.0, r0=1.0)
        @test result_gsm.status == :inner_boundary

        # Number of points might differ
        @test length(result_geo) != length(result_gsm) || result_geo.points[end] != result_gsm.points[end]
    end

    @testset "Step size effects" begin
        x, y, z = 3.0, 0.0, 0.0

        # Smaller step size should give more points
        result_fine = trace_field_line(x, y, z, t; ds=0.01, rlim=10.0, r0=1.0)
        result_coarse = trace_field_line(x, y, z, t; ds=0.1, rlim=10.0, r0=1.0)

        @test length(result_fine) > length(result_coarse)

        # Both should reach same boundary
        @test result_fine.status == result_coarse.status
    end

    @testset "FieldLineResult iteration" begin
        result = trace_field_line(3.0, 0.0, 0.0, t; rlim=10.0, r0=1.0)

        # Test iteration
        count = 0
        for point in result
            count += 1
            @test length(point) == 3
        end
        @test count == length(result)

        # Test indexing
        @test result[1] == result.points[1]
        @test result[end] == result.points[end]
    end
end

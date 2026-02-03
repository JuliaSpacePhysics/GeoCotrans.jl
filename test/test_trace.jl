using Test
using GeoCotrans
using OrdinaryDiffEq
using Dates
using LinearAlgebra

@testset "trace_field_line" begin
    t = DateTime(2020, 1, 1)

    @testset "Basic tracing from dayside equator" begin
        # Start at L=4 on dayside equator (GEO coordinates)
        pos = [4.0, 0.0, 0.0]
        sol = trace_field_line(pos, t, Tsit5(); rlim=10.0, r0=1.0)

        # Should have multiple points
        @test length(sol.u) > 10

        # First point should be the starting position
        @test sol.u[1] ≈ pos

        # Last point should be at r ≈ r0
        r_end = norm(sol.u[end])
        @test r_end ≈ 1.0 rtol=0.1

        # All intermediate radii should be >= r0 (approximately)
        @test all(u -> norm(u) >= 0.99, sol.u)
    end

    @testset "Tracing in both directions" begin
        pos = [3.0, 0.0, 0.5]

        # Trace parallel to B
        sol_fwd = trace_field_line(pos, t, Tsit5(); dir=1, rlim=10.0, r0=1.0)

        # Trace anti-parallel to B
        sol_bwd = trace_field_line(pos, t, Tsit5(); dir=-1, rlim=10.0, r0=1.0)

        # Both should reach inner boundary for dipole-like field
        r_fwd_end = norm(sol_fwd.u[end])
        r_bwd_end = norm(sol_bwd.u[end])
        @test r_fwd_end ≈ 1.0 rtol=0.1
        @test r_bwd_end ≈ 1.0 rtol=0.1

        # The endpoints should be in different hemispheres (approximately)
        z_fwd = sol_fwd.u[end][3]
        z_bwd = sol_bwd.u[end][3]
        @test sign(z_fwd) != sign(z_bwd) || abs(z_fwd) < 0.5 || abs(z_bwd) < 0.5
    end

    @testset "Outer boundary reached" begin
        # Start far from Earth with small rlim
        pos = [-5.0, 0.0, 3.0]
        sol = trace_field_line(pos, t, Tsit5(); rlim=6.0, r0=1.0, dir=1)

        # Should eventually hit some boundary
        @test length(sol.u) > 1
        r_end = norm(sol.u[end])
        # Either hit inner, outer, or lateral boundary
        @test r_end <= 6.1 || r_end >= 0.99
    end

    @testset "Vector input interface" begin
        pos = [3.0, 0.0, 0.0]
        sol = trace_field_line(pos, t, Tsit5(); rlim=10.0, r0=1.0)

        @test norm(sol.u[end]) ≈ 1.0 rtol=0.1
        @test length(sol.u) > 10
    end

    @testset "CoordinateVector input interface" begin
        pos = GEO(3.0, 0.0, 0.0)
        sol = trace_field_line(pos, t, Tsit5(); rlim=10.0, r0=1.0)

        @test norm(sol.u[end]) ≈ 1.0 rtol=0.1
        @test length(sol.u) > 10
    end

    @testset "FieldLineProblem and FieldLineCallback" begin
        pos = [3.0, 0.0, 0.0]

        # Create problem and callback separately
        prob = FieldLineProblem(pos, (0.0, 50.0), t)
        cb = FieldLineCallback(r0=1.0, rlim=10.0)

        # Solve with custom options
        sol = solve(prob, Tsit5(); callback=cb, abstol=1e-8, reltol=1e-8)

        @test length(sol.u) > 10
        @test norm(sol.u[end]) ≈ 1.0 rtol=0.1
    end

    @testset "Different solvers" begin
        pos = [3.0, 0.0, 0.0]

        # Test with different solvers
        sol_tsit5 = trace_field_line(pos, t, Tsit5(); rlim=10.0, r0=1.0)
        sol_dp5 = trace_field_line(pos, t, DP5(); rlim=10.0, r0=1.0)

        # Both should reach similar endpoint
        @test norm(sol_tsit5.u[end]) ≈ 1.0 rtol=0.1
        @test norm(sol_dp5.u[end]) ≈ 1.0 rtol=0.1
    end

    @testset "Solution iteration" begin
        pos = [3.0, 0.0, 0.0]
        sol = trace_field_line(pos, t, Tsit5(); rlim=10.0, r0=1.0)

        # Test iteration
        count = 0
        for u in sol.u
            count += 1
            @test length(u) == 3
        end
        @test count == length(sol.u)

        # Test indexing
        @test sol.u[1] == sol.u[begin]
        @test sol.u[end] == sol.u[length(sol.u)]
    end
end

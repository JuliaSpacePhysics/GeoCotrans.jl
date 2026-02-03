# Field Line Tracing

GeoCotrans.jl provides functionality to trace magnetic field lines using the IGRF model. This feature requires loading a SciML solver package (like OrdinaryDiffEq.jl).

## Basic Usage

```julia
using GeoCotrans, OrdinaryDiffEq, Dates

# Define time for field evaluation
t = DateTime(2020, 1, 1)

# Trace a field line starting from [3, 0, 0] Earth radii (GEO coordinates)
sol = trace_field_line([3.0, 0.0, 0.0], t, Tsit5(); rlim=10.0, r0=1.0)

# Access the traced points
for u in sol.u
    x, y, z = u
    println("Position: ($x, $y, $z) Re")
end
```

## API Reference

```@docs
FieldLineProblem
FieldLineCallback
trace_field_line
```

## Parameters

- `pos`: Starting position in GEO Cartesian coordinates (Earth radii)
- `t`: Time for field evaluation (DateTime, Date, or year as Float64)
- `solver`: Any SciML-compatible ODE solver (e.g., `Tsit5()`, `Vern7()`, `DP5()`)

### Keyword Arguments

- `model`: Magnetic field model (default: `IGRF()`)
- `dir`: Tracing direction, `1` for parallel to B, `-1` for anti-parallel (default: `1`)
- `r0`: Inner radial boundary in Earth radii (default: `1.0`)
- `rlim`: Outer radial boundary in Earth radii (default: `10.0`)
- `maxs`: Maximum arc length for integration (default: `100.0`)

## Plotting with Makie

Here's a complete example showing how to trace and visualize magnetic field lines using GLMakie:

```julia
using GeoCotrans, OrdinaryDiffEq, Dates
using GLMakie

# Time for field evaluation
t = DateTime(2020, 1, 1)

# Create figure
fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1],
    xlabel = "X (Re)",
    ylabel = "Y (Re)",
    zlabel = "Z (Re)",
    title = "IGRF Field Lines at $(Date(t))",
    aspect = :data
)

# Draw Earth (unit sphere)
θ_sphere = range(0, π, length=30)
φ_sphere = range(0, 2π, length=60)
x_sphere = [sin(θ) * cos(φ) for θ in θ_sphere, φ in φ_sphere]
y_sphere = [sin(θ) * sin(φ) for θ in θ_sphere, φ in φ_sphere]
z_sphere = [cos(θ) for θ in θ_sphere, φ in φ_sphere]
surface!(ax, x_sphere, y_sphere, z_sphere, color = :lightblue, alpha = 0.8)

# Trace field lines from different starting positions
L_values = [2.0, 3.0, 4.0, 5.0]  # L-shell values
colors = [:red, :orange, :green, :blue]

for (L, color) in zip(L_values, colors)
    # Trace in both directions from the equator
    for dir in [1, -1]
        sol = trace_field_line([L, 0.0, 0.0], t, Tsit5();
            dir = dir,
            r0 = 1.0,
            rlim = 10.0
        )

        # Extract coordinates
        xs = [u[1] for u in sol.u]
        ys = [u[2] for u in sol.u]
        zs = [u[3] for u in sol.u]

        # Plot the field line
        lines!(ax, xs, ys, zs, color = color, linewidth = 2,
            label = dir == 1 ? "L = $L" : nothing)
    end
end

# Add legend
axislegend(ax, position = :rt)

# Set view limits
limits!(ax, -6, 6, -6, 6, -6, 6)

fig
```

## Tracing Field Lines at Different Longitudes

```julia
using GeoCotrans, OrdinaryDiffEq, Dates
using GLMakie

t = DateTime(2020, 1, 1)

fig = Figure(size = (900, 700))
ax = Axis3(fig[1, 1],
    xlabel = "X (Re)", ylabel = "Y (Re)", zlabel = "Z (Re)",
    title = "Field Lines at Different Longitudes",
    aspect = :data
)

# Draw Earth
θ_sphere = range(0, π, length=30)
φ_sphere = range(0, 2π, length=60)
x_sphere = [sin(θ) * cos(φ) for θ in θ_sphere, φ in φ_sphere]
y_sphere = [sin(θ) * sin(φ) for θ in θ_sphere, φ in φ_sphere]
z_sphere = [cos(θ) for θ in θ_sphere, φ in φ_sphere]
surface!(ax, x_sphere, y_sphere, z_sphere, color = :lightblue, alpha = 0.6)

# Trace field lines at different longitudes
L = 4.0
longitudes = range(0, 2π, length=9)[1:8]  # 8 longitudes

for φ in longitudes
    # Starting position on the magnetic equator
    x0 = L * cos(φ)
    y0 = L * sin(φ)
    z0 = 0.0

    for dir in [1, -1]
        sol = trace_field_line([x0, y0, z0], t, Tsit5();
            dir = dir, r0 = 1.0, rlim = 10.0)

        xs = [u[1] for u in sol.u]
        ys = [u[2] for u in sol.u]
        zs = [u[3] for u in sol.u]

        lines!(ax, xs, ys, zs, color = :navy, linewidth = 1.5)
    end
end

limits!(ax, -5, 5, -5, 5, -5, 5)
fig
```

## Advanced: Using FieldLineProblem Directly

For more control over the integration, you can use the lower-level API:

```julia
using GeoCotrans, OrdinaryDiffEq, Dates

t = DateTime(2020, 1, 1)
pos = [4.0, 0.0, 0.0]

# Create the ODE problem
prob = FieldLineProblem(pos, (0.0, 100.0), t; dir = 1)

# Create custom callback with different boundaries
cb = FieldLineCallback(r0 = 1.0, rlim = 15.0, xlim = 25.0)

# Solve with custom tolerances and dense output
sol = solve(prob, Vern7();
    callback = cb,
    abstol = 1e-10,
    reltol = 1e-10,
    dense = true
)

# With dense output, you can interpolate at any arc length
s_values = range(0, sol.t[end], length=100)
positions = sol.(s_values)
```

## Dipole Field Comparison

You can visualize how IGRF differs from a simple dipole by tracing field lines:

```julia
using GeoCotrans, OrdinaryDiffEq, Dates
using GLMakie

t = DateTime(2020, 1, 1)

fig = Figure(size = (1000, 500))

# Left panel: Meridional plane (X-Z)
ax1 = Axis(fig[1, 1],
    xlabel = "X (Re)", ylabel = "Z (Re)",
    title = "Field Lines in Meridional Plane",
    aspect = DataAspect()
)

# Draw Earth cross-section
θ_circle = range(0, 2π, length=100)
lines!(ax1, cos.(θ_circle), sin.(θ_circle), color = :lightblue, linewidth = 3)

# Trace field lines
for L in [2.0, 3.0, 4.0, 5.0, 6.0]
    for dir in [1, -1]
        sol = trace_field_line([L, 0.0, 0.0], t, Tsit5();
            dir = dir, r0 = 1.0, rlim = 10.0)

        xs = [u[1] for u in sol.u]
        zs = [u[3] for u in sol.u]

        lines!(ax1, xs, zs, color = :darkblue, linewidth = 1.5)
    end
end

limits!(ax1, -7, 7, -7, 7)

# Right panel: Equatorial plane (X-Y)
ax2 = Axis(fig[1, 2],
    xlabel = "X (Re)", ylabel = "Y (Re)",
    title = "Field Lines in Equatorial Plane",
    aspect = DataAspect()
)

# Draw Earth cross-section
lines!(ax2, cos.(θ_circle), sin.(θ_circle), color = :lightblue, linewidth = 3)

# Starting points slightly off equator to see field line structure
for L in [3.0, 4.0, 5.0]
    for φ in range(0, 2π, length=13)[1:12]
        x0 = L * cos(φ)
        y0 = L * sin(φ)

        # Start slightly above equator
        sol = trace_field_line([x0, y0, 0.1], t, Tsit5();
            dir = 1, r0 = 1.0, rlim = 10.0)

        xs = [u[1] for u in sol.u]
        ys = [u[2] for u in sol.u]

        lines!(ax2, xs, ys, color = :darkgreen, linewidth = 1)
    end
end

limits!(ax2, -7, 7, -7, 7)

fig
```

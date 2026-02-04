# Field Line Tracing

```@docs; canonical=false
GeoCotrans.FieldLineTracing
```

## Basic Usage

```julia
using GeoCotrans, OrdinaryDiffEqTsit5, Dates

# Define time for field evaluation
t = DateTime(2020, 1, 1)

# Trace a field line starting from [3, 0, 0] Earth radii (GEO coordinates)
sol = trace(GEO(3.0, 0.0, 0.0), t, Tsit5())
```

## API Reference

```@docs; canonical=false
FieldLineProblem
FieldLineCallback
trace
```

## Plotting with Makie

Here's a complete example showing how to trace and visualize magnetic field lines using GLMakie:

```julia
using GeoCotrans, OrdinaryDiffEqTsit5, Dates
using CairoMakie

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
        sol = trace(GEO(L, 0.0, 0.0), t, Tsit5(); dir = dir)
        # Plot the field line
        plot!(ax, sol; color = color, linewidth = 2, idxs = (1, 2, 3),
            label = dir == 1 ? "L = $L" : nothing)
    end
end

axislegend(ax, position = :rt) # Add legend
limits!(ax, -6, 6, -6, 6, -6, 6) # Set view limits

fig
```

## 2D Visualization

```julia
using GeoCotrans, OrdinaryDiffEqTsit5, Dates
using CairoMakie

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
        sol = trace(GEO(L, 0.0, 0.0), t, Tsit5(); dir = dir)
        lines!(ax1, sol, idxs = (1, 3), color = :darkblue, linewidth = 1.5)
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
        y0, x0 = L .* sincos(φ)
        # Start slightly above equator
        sol = trace(GEO(x0, y0, 0.1), t, Tsit5(); dir = 1)
        lines!(ax2, sol, idxs = (1, 2), color = :darkgreen, linewidth = 1)
    end
end

limits!(ax2, -7, 7, -7, 7)

fig
```

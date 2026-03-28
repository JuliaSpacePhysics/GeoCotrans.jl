# GeoCotrans.jl Development Workflow

Guide for common development tasks in this Julia package for geocentric coordinate transformations and IGRF field modeling.

## Setup

Instantiate dependencies:
```
julia --project -e "using Pkg; Pkg.instantiate()"
```

## Running Tests

```
julia --project --threads=auto test/runtests.jl
```

Tests use `TestItems.jl` / `TestItemRunner`. On stable Julia, JET static analysis also runs. To run a specific test item, use `@run_package_tests filter=ti->contains(ti.name, "keyword")`.

## Building Docs

```
cd docs && julia --project make.jl
```

Docs use Documenter.jl with doctests enabled. The doc site covers coordinate system descriptions, validation comparisons (vs IRBEM/PySPEDAS), and field line tracing examples.

## Key Source Layout

- `src/cotrans/` — coordinate transformation pairs (GEI, GEO, GSM, GSE, MAG, SM)
- `src/igrf.jl` + `src/spherical_harmonics.jl` — IGRF-14 field model
- `src/types.jl` — `CoordinateVector{Frame, Representation}` type
- `src/trace.jl` + `ext/GeoCotransSciMLExt.jl` — field line tracing via ODE
- `ext/GeoCotransDimensionalDataExt.jl` — DimensionalData integration

## Git Workflow

Development branch: `claude/repo-workflow-skill-CwNSO`

```bash
# Stage and commit
git add <files>
git commit -m "description"

# Push
git push -u origin claude/repo-workflow-skill-CwNSO
```

CI runs on push/PR: tests on Julia stable/LTS/pre-release, coverage to Codecov, docs deploy on main.

## Adding a Transformation

1. Add `src/cotrans/<from>2<to>.jl` implementing the rotation matrix
2. Register the pair in `src/cotrans/cotrans.jl` via the macro system
3. Export the new function from `src/GeoCotrans.jl`
4. Add test items in `test/`
5. Add validation comparison in `docs/src/coords.md`

# Newton Fractal

An interactive **Newton fractal** explorer written in pure Julia. For every starting
point `z = a + bi` in a region of the complex plane, it runs Newton's method on a
chosen polynomial and colors the pixel by **which root the iteration converges to**
(its *basin of attraction*). Shading by the number of iterations reveals the
infinitely detailed boundary between basins.

The math: by the Fundamental Theorem of Algebra a degree-`n` polynomial has `n`
roots in `ℂ`, and Newton's iteration

```
z_{k+1} = z_k − p(z_k) / p'(z_k)
```

extends to the complex plane unchanged. Each pixel is one starting value; its final
landing point is matched to the nearest root.

## Requirements

- [Julia](https://julialang.org/downloads/) **1.9 or newer**
- A desktop with a GPU/OpenGL (GLMakie opens a native window)
- Dependencies (`GLMakie`, `Polynomials`) install automatically on first run

## Quick start

From this folder:

```bash
julia run.jl
```

That's the whole thing. On the **first** run `run.jl` will:

1. activate the local project (`Project.toml`)
2. download and install `GLMakie` + `Polynomials`
3. relaunch itself with `--threads=auto` so the grid is computed in parallel
4. open the interactive window

The first launch is slow (Julia precompiles GLMakie — a few minutes). Every run
after that opens in seconds.

> Prefer the REPL? Do `using Pkg; Pkg.activate("."); Pkg.instantiate()` once, then
> `include("newton_fractal.jl")`.

## Using the window

| Control | What it does |
|---|---|
| **Polinômio** | Pick a preset, defined either by its roots or its coefficients |
| **Paleta** | Color map for the basins |
| **Máx. iterações (M)** | Iteration cap per pixel (higher = sharper boundaries, slower) |
| **Resolução** | Grid size (higher = crisper, slower) |
| **sombrear** | Toggle iteration-count shading (the fractal "glow") |
| **Renderizar região (zoom)** | Scroll/drag to zoom the axis, then click this to **recompute** the fractal at that zoom level — true deep zoom, not a stretched image |
| **Reset zoom** | Back to the default `[-2, 2] × [-2, 2]` view |
| **Raízes (text field)** | Type roots, comma-separated, and press Enter — e.g. `1, -0.5+0.87im, -0.5-0.87im` |
| **Coefs (text field)** | Type coefficients `a₀, a₁, …, aₙ` (lowest degree first) and press Enter — e.g. `-1, 0, 0, 1` for `z³ − 1` |

The imaginary unit is Julia's `im` (so `2im`, not `2i`). Format errors show up in the info label instead of crashing.

Things to look for: in `z³ − 2z + 2` the black patches belong to **no** basin —
those starting points fall into an attracting cycle and never converge (Smale's
classic example). And on any boundary between two colors, zooming in always reveals
the third color squeezed in between (the Wada property) — that's what makes the
set a fractal.

## Project structure

```
fractals/
├── run.jl              # launcher: installs deps, sets threads, opens the window
├── newton_fractal.jl   # the app (numerics + GLMakie UI)
├── Project.toml        # declared dependencies
└── README.md
```

## Customizing

- **Add a polynomial:** append to the `PRESETS` list in `newton_fractal.jl`, either
  as `:roots => [r1, r2, ...]` or `:coeffs => [a0, a1, ..., an]` (lowest degree first) —
  or just type roots/coefficients directly in the window's text fields (press Enter).
- **Add a palette:** add any Makie colormap name (a `Symbol`) to `PALETTES`.

## Troubleshooting

- **First run takes minutes / seems stuck:** that's GLMakie precompiling. Let it finish; subsequent runs are fast.
- **Rendering feels slow:** make sure threads kicked in — the console prints
  `Threads in use: N`. If it says `1`, run `julia --threads=auto run.jl` directly.
- **No window appears (headless/SSH):** GLMakie needs a real display. On a remote
  machine, swap `using GLMakie` for `using CairoMakie` to render to a file instead.

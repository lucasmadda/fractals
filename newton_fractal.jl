# ======================================================================
#  Fractal de Newton — interativo, Julia puro (GLMakie + Polynomials)
#  repo: fractals          autor: Lucas
# ----------------------------------------------------------------------
#  Easiest way to run:   julia run.jl     (installs deps + threads + opens)
#
#  From the REPL instead:
#     julia> using Pkg; Pkg.activate("."); Pkg.instantiate()
#     julia> include("newton_fractal.jl")
#
#  Window features:
#    • polynomial by ROOTS or COEFFICIENTS — presets OR typed directly
#      in the text fields (comma-separated, imaginary unit `im`, Enter)
#    • color palettes + shading by iteration count (the fractal "glow")
#    • zoom: scroll/drag the axis, then click "Renderizar região"
#      to RE-COMPUTE the fractal in that window (true deep zoom)
# ======================================================================

using GLMakie, Polynomials

# ----------------------------------------------------------------------
# 1) Numerical core: iterate Newton over a grid of the complex plane.
#    Returns the basin (root index) and speed (iters/M) per pixel.
# ----------------------------------------------------------------------
function newton_basins(P::Polynomial, dP::Polynomial, rs::Vector{<:Complex};
                       xlims, ylims, res::Int, maxiter::Int, tol::Float64)
    xs = range(xlims[1], xlims[2]; length = res)
    ys = range(ylims[1], ylims[2]; length = res)
    basin = Matrix{Int}(undef, res, res)
    speed = Matrix{Float64}(undef, res, res)

    Threads.@threads for j in 1:res
        @inbounds for i in 1:res
            z = complex(xs[i], ys[j])
            k = 0
            while k < maxiter
                d = dP(z)
                abs(d) < 1e-14 && break          # derivative ~ 0: critical point
                znew = z - P(z) / d
                k += 1
                if abs(znew - z) < tol           # converged
                    z = znew
                    break
                end
                z = znew
            end
            best, bd = 1, Inf                    # argmin |z - rj|
            for (m, r) in enumerate(rs)
                dd = abs2(z - r)
                if dd < bd
                    bd, best = dd, m
                end
            end
            basin[i, j] = best
            speed[i, j] = k / maxiter
        end
    end
    return basin, speed
end

# ----------------------------------------------------------------------
# 2) basins + speed  ->  RGB image
#    base color per root (sampled from palette) x glow(k) = (1 - k/M)^gamma
#    NOTE: we read the .r/.g/.b fields directly — Makie does NOT export
#    the Colors accessors red()/green()/blue(), so calling them errors.
# ----------------------------------------------------------------------
function basins_to_image(basin, speed, n::Int, palette::Symbol; shade::Bool, γ = 0.45)
    grad = cgrad(palette)
    base = [grad[clamp((m - 1) / max(n - 1, 1), 0, 1)] for m in 1:n]
    r1, r2 = size(basin)
    img = Matrix{RGBf}(undef, r1, r2)
    @inbounds for j in 1:r2, i in 1:r1
        c = base[basin[i, j]]
        f = shade ? (1 - speed[i, j])^γ : 1.0
        img[i, j] = RGBf(c.r * f, c.g * f, c.b * f)
    end
    return img
end

# ----------------------------------------------------------------------
# 3) Polynomial builders + input parsing
# ----------------------------------------------------------------------
poly_from_roots(rs)  = fromroots(ComplexF64.(rs))
function poly_from_coeffs(c)                          # c = [a0, a1, ..., an]
    P = Polynomial(ComplexF64.(c))
    return P, ComplexF64.(roots(P))
end

# "1, -0.5+0.87im, 2im"  ->  Vector{ComplexF64}   (imaginary unit: im)
parse_clist(s::AbstractString) =
    [parse(ComplexF64, replace(t, " " => "")) for t in split(s, ",") if !isempty(strip(t))]

# ----------------------------------------------------------------------
# 4) Presets  (:roots => root list  |  :coeffs => [a0..an])
# ----------------------------------------------------------------------
const PRESETS = [
    "z³ − 1"                 => (:roots,  [cis(2π*k/3) for k in 0:2]),
    "z⁴ − 1"                 => (:roots,  [cis(2π*k/4) for k in 0:3]),
    "z⁵ − 1"                 => (:roots,  [cis(2π*k/5) for k in 0:4]),
    "z⁶ − 1"                 => (:roots,  [cis(2π*k/6) for k in 0:5]),
    "z³ − 2z + 2 (chaotic)"  => (:coeffs, [2.0, -2.0, 0.0, 1.0]),
    "z⁶ + z³ − 1"            => (:coeffs, [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]),
    "z⁸ + 15z⁴ − 16"         => (:coeffs, [-16.0, 0,0,0, 15.0, 0,0,0, 1.0]),
]
const PALETTES = [:viridis, :plasma, :turbo, :magma, :inferno, :rainbow]
const LOOKUP   = Dict(PRESETS)

# ----------------------------------------------------------------------
# 5) Interactive app  (returns the Figure; does NOT display it)
# ----------------------------------------------------------------------
function launch()
    fig = Figure(size = (1200, 840))
    ax  = Axis(fig[1, 2]; aspect = DataAspect(),
               title = "Fractal de Newton", xlabel = "Re", ylabel = "Im")

    ctrl = fig[1, 1] = GridLayout(tellheight = false, valign = :top)
    Label(ctrl[1, 1:2], "Controles"; font = :bold, fontsize = 18)
    Label(ctrl[2, 1:2], "Polinômio");          menu_preset = Menu(ctrl[3, 1:2]; options = first.(PRESETS))
    Label(ctrl[4, 1:2], "Paleta");             menu_pal    = Menu(ctrl[5, 1:2]; options = string.(PALETTES))
    Label(ctrl[6, 1:2], "Máx. iterações (M)"); sl_iter     = Slider(ctrl[7, 1:2]; range = 20:10:300, startvalue = 70)
    Label(ctrl[8, 1:2], "Resolução");          sl_res      = Slider(ctrl[9, 1:2]; range = 200:100:1400, startvalue = 600)
    tg_shade = Toggle(ctrl[10, 1]; active = true); Label(ctrl[10, 2], "sombrear (brilho)")
    btn_render = Button(ctrl[11, 1:2]; label = "Renderizar região (zoom)")
    btn_reset  = Button(ctrl[12, 1:2]; label = "Reset zoom")
    info       = Label(ctrl[13, 1:2], "—"; tellwidth = false)
    Label(ctrl[14, 1:2], "Raízes (Enter aplica):")
    tb_roots  = Textbox(ctrl[15, 1:2]; placeholder = "1, -0.5+0.87im, -0.5-0.87im", width = 230)
    Label(ctrl[16, 1:2], "Coefs a₀,…,aₙ (Enter):")
    tb_coeffs = Textbox(ctrl[17, 1:2]; placeholder = "-1, 0, 0, 1", width = 230)

    state   = Ref{Any}(nothing)                 # (P, dP, rs)
    DEFAULT = ((-2.0, 2.0), (-2.0, 2.0))
    cur     = Ref(DEFAULT)

    function set_poly!(name)
        kind, data = LOOKUP[name]
        if kind === :roots
            P, rs = poly_from_roots(data), ComplexF64.(data)
        else
            P, rs = poly_from_coeffs(data)
        end
        state[] = (P, derivative(P), rs)
        info.text = "grau $(degree(P)) · $(length(rs)) raízes"
    end

    function render!(xlims, ylims)
        state[] === nothing && return
        P, dP, rs = state[]
        basin, speed = newton_basins(P, dP, rs;
            xlims = xlims, ylims = ylims, res = Int(sl_res.value[]),
            maxiter = Int(sl_iter.value[]), tol = 1e-6)
        img = basins_to_image(basin, speed, length(rs),
            Symbol(menu_pal.selection[]); shade = tg_shade.active[])
        empty!(ax)
        image!(ax, xlims[1] .. xlims[2], ylims[1] .. ylims[2], img)
        limits!(ax, xlims[1], xlims[2], ylims[1], ylims[2])
        cur[] = (xlims, ylims)
    end

    # custom polynomial from the text fields ----------------------------
    function set_custom!(P, rs, tag)
        if degree(P) < 1 || isempty(rs)
            info.text = "preciso de grau ≥ 1"
            return
        end
        state[] = (P, derivative(P), rs)
        info.text = "grau $(degree(P)) · $(length(rs)) raízes · $tag"
        render!(DEFAULT...)
    end
    on(tb_roots.stored_string) do s
        s === nothing && return
        try
            rs = parse_clist(s)
            set_custom!(poly_from_roots(rs), ComplexF64.(rs), "raízes")
        catch
            info.text = "raízes: formato 1, -0.5+0.87im, … (use im)"
        end
    end
    on(tb_coeffs.stored_string) do s
        s === nothing && return
        try
            P, rs = poly_from_coeffs(parse_clist(s))
            set_custom!(P, rs, "coefs")
        catch
            info.text = "coefs: formato a₀, a₁, … (use im)"
        end
    end

    # other events -------------------------------------------------------
    on(menu_preset.selection) do name
        name === nothing && return
        set_poly!(name); render!(DEFAULT...)
    end
    for o in (menu_pal.selection, sl_iter.value, sl_res.value, tg_shade.active)
        on(o) do _
            state[] === nothing && return
            render!(cur[]...)
        end
    end
    on(btn_render.clicks) do _                  # recompute in the zoomed window
        r = ax.finallimits[]
        render!((r.origin[1], r.origin[1] + r.widths[1]),
                (r.origin[2], r.origin[2] + r.widths[2]))
    end
    on(btn_reset.clicks) do _; render!(DEFAULT...); end

    set_poly!(first(PRESETS)[1]); render!(DEFAULT...)
    return fig
end

# ----------------------------------------------------------------------
# 6) Show it. In the REPL: opens and returns. As a script (run.jl):
#    opens and blocks until the window is closed.
# ----------------------------------------------------------------------
let screen = display(launch())
    isinteractive() || wait(screen)
end

# ----------------------------------------------------------------------
# NOTE — scroll-zoom that recomputes automatically: call render!(cur[]...)
# inside  on(events(ax.scene).scroll)  with the axis' current limits.
# ----------------------------------------------------------------------

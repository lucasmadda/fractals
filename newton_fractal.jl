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
#    • UI language switch (Português / English) at the top
#    • polynomial by ROOTS or COEFFICIENTS — presets OR typed directly
#      (comma-separated, imaginary unit `im`, press Enter)
#    • roots drawn as white dots on top of their basins
#    • live slider readouts, palettes, iteration shading
#    • zoom: scroll/drag the axis, then "Renderizar região" RE-COMPUTES
#      the fractal in that window (true deep zoom)
# ======================================================================

using GLMakie, Polynomials

GLMakie.activate!(title = "Newton fractal")

# ----------------------------------------------------------------------
# i18n — every UI string in both languages
# ----------------------------------------------------------------------
const STR = Dict(
    :pt => Dict(
        :controls => "Controles",
        :poly     => "Polinômio",
        :pal      => "Paleta",
        :iter     => "Máx. iterações (M)",
        :res      => "Resolução",
        :shade    => "sombrear (brilho)",
        :render   => "Renderizar região (zoom)",
        :reset    => "Reset zoom",
        :roots_in => "Raízes (Enter aplica):",
        :coefs_in => "Coefs a₀,…,aₙ (Enter):",
        :title    => "Fractal de Newton",
        :degree   => "grau",
        :roots    => "raízes",
        :custom   => "personalizado",
        :needdeg  => "preciso de grau ≥ 1",
        :err_r    => "raízes: formato 1, -0.5+0.87im, … (use im)",
        :err_c    => "coefs: formato a₀, a₁, … (use im)",
    ),
    :en => Dict(
        :controls => "Controls",
        :poly     => "Polynomial",
        :pal      => "Palette",
        :iter     => "Max iterations (M)",
        :res      => "Resolution",
        :shade    => "shading (glow)",
        :render   => "Render region (zoom)",
        :reset    => "Reset zoom",
        :roots_in => "Roots (Enter applies):",
        :coefs_in => "Coeffs a₀,…,aₙ (Enter):",
        :title    => "Newton fractal",
        :degree   => "degree",
        :roots    => "roots",
        :custom   => "custom",
        :needdeg  => "need degree ≥ 1",
        :err_r    => "roots: format 1, -0.5+0.87im, … (use im)",
        :err_c    => "coeffs: format a₀, a₁, … (use im)",
    ),
)
const LANGS = ["Português" => :pt, "English" => :en]

# ----------------------------------------------------------------------
# 1) Numerical core: iterate Newton over a grid of the complex plane.
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
#    (.r/.g/.b field access — Makie does not export red()/green()/blue())
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
    "z³ − 2z + 2"            => (:coeffs, [2.0, -2.0, 0.0, 1.0]),
    "z⁶ + z³ − 1"            => (:coeffs, [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0]),
    "z⁸ + 15z⁴ − 16"         => (:coeffs, [-16.0, 0,0,0, 15.0, 0,0,0, 1.0]),
]
const PALETTES = [:viridis, :plasma, :turbo, :magma, :inferno, :rainbow]
const LOOKUP   = Dict(PRESETS)

# ----------------------------------------------------------------------
# 5) Interactive app  (returns the Figure; does NOT display it)
# ----------------------------------------------------------------------
function launch()
    lang = Ref(:pt)
    t(k) = STR[lang[]][k]

    fig = Figure(size = (1240, 860))
    ax  = Axis(fig[1, 2]; aspect = DataAspect(), xlabel = "Re", ylabel = "Im")

    ctrl = fig[1, 1] = GridLayout(tellheight = false, valign = :top)
    colsize!(fig.layout, 1, Fixed(300))      # fixed panel width: nothing overflows

    lbl_head   = Label(ctrl[1, 1:2], ""; font = :bold, fontsize = 18, halign = :left)
    menu_lang  = Menu(ctrl[2, 1:2]; options = first.(LANGS))
    lbl_poly   = Label(ctrl[3, 1:2], ""; halign = :left)
    menu_preset = Menu(ctrl[4, 1:2]; options = first.(PRESETS))
    lbl_pal    = Label(ctrl[5, 1:2], ""; halign = :left)
    menu_pal   = Menu(ctrl[6, 1:2]; options = string.(PALETTES))
    lbl_iter   = Label(ctrl[7, 1], ""; halign = :left)
    sl_iter    = Slider(ctrl[8, 1:2]; range = 20:10:300, startvalue = 70)
    Label(ctrl[7, 2], lift(string, sl_iter.value); halign = :right)   # live readout
    lbl_res    = Label(ctrl[9, 1], ""; halign = :left)
    sl_res     = Slider(ctrl[10, 1:2]; range = 200:100:1400, startvalue = 600)
    Label(ctrl[9, 2], lift(string, sl_res.value); halign = :right)    # live readout
    tg_shade   = Toggle(ctrl[11, 1]; active = true, halign = :left)
    lbl_shade  = Label(ctrl[11, 2], ""; halign = :left)
    btn_render = Button(ctrl[12, 1:2]; label = " ")
    btn_reset  = Button(ctrl[13, 1:2]; label = " ")
    info       = Label(ctrl[14, 1:2], "—"; halign = :left, tellwidth = false)
    lbl_roots  = Label(ctrl[15, 1:2], ""; halign = :left)
    tb_roots   = Textbox(ctrl[16, 1:2]; placeholder = "1, 2im, -1-1im", width = 280)
    lbl_coefs  = Label(ctrl[17, 1:2], ""; halign = :left)
    tb_coeffs  = Textbox(ctrl[18, 1:2]; placeholder = "-1, 0, 0, 1", width = 280)
    rowgap!(ctrl, 10)

    state    = Ref{Any}(nothing)                # (P, dP, rs)
    DEFAULT  = ((-2.0, 2.0), (-2.0, 2.0))
    cur      = Ref(DEFAULT)
    polyname = Ref(first(PRESETS)[1])
    inforef  = Ref((deg = 0, n = 0, custom = false))

    function refresh_text!()                    # re-applies every string in lang[]
        lbl_head.text  = t(:controls)
        lbl_poly.text  = t(:poly)
        lbl_pal.text   = t(:pal)
        lbl_iter.text  = t(:iter)
        lbl_res.text   = t(:res)
        lbl_shade.text = t(:shade)
        btn_render.label = t(:render)
        btn_reset.label  = t(:reset)
        lbl_roots.text = t(:roots_in)
        lbl_coefs.text = t(:coefs_in)
        s = inforef[]
        ax.title  = t(:title) * " — " * (s.custom ? t(:custom) : polyname[])
        info.text = s.deg == 0 ? "—" :
            "$(t(:degree)) $(s.deg) · $(s.n) $(t(:roots))" * (s.custom ? " · $(t(:custom))" : "")
    end

    function set_poly!(name)
        kind, data = LOOKUP[name]
        P, rs = kind === :roots ? (poly_from_roots(data), ComplexF64.(data)) :
                                  poly_from_coeffs(data)
        state[]    = (P, derivative(P), rs)
        polyname[] = name
        inforef[]  = (deg = degree(P), n = length(rs), custom = false)
        refresh_text!()
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
        scatter!(ax, real.(rs), imag.(rs); color = :white,        # the stars of the show
                 strokecolor = :black, strokewidth = 1.5, markersize = 11)
        limits!(ax, xlims[1], xlims[2], ylims[1], ylims[2])
        cur[] = (xlims, ylims)
    end

    # custom polynomial from the text fields ----------------------------
    function set_custom!(P, rs)
        if degree(P) < 1 || isempty(rs)
            info.text = t(:needdeg)
            return
        end
        state[]   = (P, derivative(P), rs)
        inforef[] = (deg = degree(P), n = length(rs), custom = true)
        refresh_text!()
        render!(DEFAULT...)
    end
    on(tb_roots.stored_string) do s
        s === nothing && return
        try
            rs = parse_clist(s)
            set_custom!(poly_from_roots(rs), ComplexF64.(rs))
        catch
            info.text = t(:err_r)
        end
    end
    on(tb_coeffs.stored_string) do s
        s === nothing && return
        try
            P, rs = poly_from_coeffs(parse_clist(s))
            set_custom!(P, rs)
        catch
            info.text = t(:err_c)
        end
    end

    # other events -------------------------------------------------------
    on(menu_lang.selection) do sel
        sel === nothing && return
        lang[] = Dict(LANGS)[sel]
        refresh_text!()
    end
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

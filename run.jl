#!/usr/bin/env julia
# ======================================================================
#  run.jl — one-command launcher for the Newton fractal app.
#
#  Just run:   julia run.jl
#
#  It will:
#    1. activate the local project (Project.toml in this folder)
#    2. install GLMakie + Polynomials the first time (and only then)
#    3. restart itself with --threads=auto so the grid renders in parallel
#       (thread count can't be changed after Julia starts, hence the restart)
#    4. open the interactive window
# ======================================================================

import Pkg
const HERE = @__DIR__
Pkg.activate(HERE)

# --- 1 & 2) make sure dependencies are present ------------------------
try
    Pkg.instantiate()                       # installs whatever is missing
catch err
    @warn "instantiate failed, adding deps explicitly" err
    Pkg.add(["GLMakie", "Polynomials"])
    Pkg.instantiate()
end

# --- 3) relaunch with all CPU threads if started single-threaded ------
if Threads.nthreads() == 1 && get(ENV, "NEWTON_RELAUNCHED", "") == ""
    @info "Restarting with --threads=auto for parallel rendering…"
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    cmd = addenv(`$julia --threads=auto --startup-file=no --project=$HERE $(@__FILE__)`,
                 "NEWTON_RELAUNCHED" => "1")
    ok = success(cmd)   # child's errors print once; no duplicate stacktrace here
    exit(ok ? 0 : 1)
end

# --- 4) launch -------------------------------------------------------
@info "Threads in use: $(Threads.nthreads())"
include(joinpath(HERE, "newton_fractal.jl"))

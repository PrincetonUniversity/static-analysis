# StaticAnalysis
# ==============


# Utils
# -----

"Vec2 is a 2D vector."
immutable Vec2
    x::Float64
    y::Float64
end

"A State contains the position and velocity of a particle."
immutable State
    pos::Vec2
    vel::Vec2
end

import Base: +, -, *, /
for op in [:+, :-, :*, :/]
    @eval $op(u::Vec2, v::Vec2) = Vec2($op(u.x, v.x), $op(u.y, v.y))
    @eval $op(x::Real, u::Vec2) = Vec2($op(x, u.x), $op(x, u.y))
    @eval $op(u::Vec2, x::Real) = Vec2($op(u.x, x), $op(u.y, x))
end

import Base: norm, angle, dot, cross, mod
norm(u::Vec2) = hypot(u.x, u.y)
dist(u::Vec2, v::Vec2) = norm(u - v)
unit(u::Vec2) = u / norm(u)
angle(u::Vec2) = atan2(u.y, u.x)
angle(u::Vec2, v::Vec2) = mod2pi(angle(u) - angle(v) + 3π) - π
dot(u::Vec2, v::Vec2) = u.x * v.x + u.y * v.y
cross(u::Vec2, v::Vec2) = u.x * v.y - u.y * v.x
mod(u::Vec2, m) = Vec2(mod(u.x, m), mod(u.y, m))

rotate(u::Vec2, θ::Real) = Vec2(u.x * cos(θ) - u.y * sin(θ), u.x * sin(θ) + u.y * cos(θ))


import Base: mean
function mean(p::Vector{State})
    n = length(p)
    mx, my, vx, vy = 0.0, 0.0, 0.0, 0.0
    @inbounds for i in eachindex(p)
        mx += p[i].pos.x
        my += p[i].pos.y
        vx += p[i].vel.x
        vy += p[i].vel.y
    end
    State(Vec2(mx / n, my / n), Vec2(vx / n, vy / n))
end


# Packages
# --------

# Data
using DataFrames

# Graphics
using Gadfly
using Compose
using Colors

# I/O
using HDF5
using MAT
using JLD
using Formatting

# Graphs
using LightGraphs


# Projects
# --------
#
# Each file contains *step* functions to generate or load data,
# compute derived data (typically a DataFrame), and generate plots.

"""A `Project` derives data using *step* functions.

It is initialized with a specific configuration and
also contains dependency information for all steps.
"""
type Project
    uid::AbstractString
    conf::Dict{Symbol, Any}
    step::Dict{Function, Any}
    deps::Dict{Function, Set{Function}}
    Project(uid; config...) = new(uid, Dict{Symbol, Any}(config), Dict(), Dict())
    Project(p::Project, uid; config...) = new(uid, merge(p.conf, Dict{Symbol, Any}(config)), Dict(), Dict())
end

"""Reset a step in a `Project` as well as all steps downstream."""
function reset!(p::Project, step)
    dir = joinpath(data_path, "Projects", p.uid)
    bak = joinpath(dir, "backup")
    mkpath(bak)
    try
        mv(joinpath(dir, "$(step).jld"),
           joinpath(bak, "$(step).jld"),
           remove_destination=true)
    end
    delete!(p.step, step)
    reset_from!(p, step)
end

"""Reset all steps downstream of a step in a `Project`."""
function reset_from!(p::Project, step)
    for (s, deps) in p.deps
        step in deps && reset!(p, s)
    end
end

"""Rebuild a `Project` with missing steps."""
function rebuild!(p::Project)
    for (step, _) in p.deps
        get!(()->step(p), p.step, step)
    end
    println("Done")
end

"""`@step` declares a *step* function, which takes a `Project` and
populates the corresponding entry.

The function may only use the entries of the `Project` indicated as
dependencies in an optional vector of step functions immediatly
following the `@step` declaration.

### Example

    @step function step_one(p::Project)
        true
    end

    @step function step_two(p::Project)
        12
    end

    @step [step_one, step_two] function step_three(p::Project)
        p[step_one] ? 42 : p[step_two] + 1
    end
"""
macro step(args...)
    if length(args) == 1
        deps, ex = Symbol[], args[1]
    elseif length(args) == 2
        deps, ex = args
    else
        error("Usage: @step [deps] func")
    end
    name, p = ex.args[1].args[1:2]
    ex.args[1].args[1] = esc(ex.args[1].args[1])
    ex.args[2] = quote
        dir = joinpath(data_path, "Projects", $p.uid)
        get!($p.deps, $name, Set{Function}())
        for dep in $deps
            push!($p.deps[$name], dep)
            get!($p.step, dep) do
                data = dep($p)
                JLD.save(joinpath(dir, "$(dep).jld"), "data", data)
                data
            end
        end
        file = joinpath(dir, "$($name).jld")
        if isfile(file)
            println("Loading ", $name, ":")
            $p.step[$name] = JLD.load(file, "data")
        else
            mkpath(dir)
            println("Running ", $name, ":")
            $p.step[$name] = begin $(ex.args[2]) end
            JLD.save(file, "data", $p.step[$name])
            $p.step[$name]
        end
    end
    ex
end


# Colormaps and othe miscellaneous things.
include("misc.jl")

# Configuration(s).
include("config.jl")

# Consolidate and reformat datasets for ellipswarm.
include("io.jl")

# Compute α-shapes and distances from the edge of the α-shapes.
include("alpha.jl")

# Compute personal and social information from ellipswarm/static.
include("info.jl")

# Analyze detection count in and around the schools.
include("detections.jl")

# Strip detections based on visual zones (binocular, blind…)
include("zones.jl")

# Test cut transversal connections
include("cut.jl")

# Analyze structure of visual interaction networks.
include("structure.jl")

# Analyze correlation of information in detections.
include("spread.jl")

# Analyze effect of changing interindividual distances.
include("rescaling.jl")

# Analyze effect of turbidity on schooling behavior.
include("turbidity.jl")

# Analyze effect of predator cooperation on risk.
include("predator.jl")


nothing

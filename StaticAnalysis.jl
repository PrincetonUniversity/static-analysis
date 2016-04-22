# StaticAnalysis
# ==============


# Settings
# --------

data_path = joinpath(homedir(), "data", "static-analysis")
plot_path = joinpath(homedir(), "plots", "static-analysis")

α_radius = -0.2

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

import Base: norm, angle, dot, cross
norm(u::Vec2) = hypot(u.x, u.y)
dist(u::Vec2, v::Vec2) = norm(u - v)
unit(u::Vec2) = u / norm(u)
angle(u::Vec2) = atan2(u.y, u.x)
angle(u::Vec2, v::Vec2) = mod2pi(angle(u) - angle(v) + 3π) - π
dot(u::Vec2, v::Vec2) = u.x * v.x + u.y * v.y
cross(u::Vec2, v::Vec2) = u.x * v.y - u.y * v.x

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

# Graphs
using LightGraphs


# Projects
# --------
#
# Each file contains functions to generate or load data,
# compute derived data (typically a DataFrame), and generate plots.

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

# Analyse structure of visual interaction networks.
include("structure.jl")

# Analyse correlation of information in detections.
include("spread.jl")

# Analyse effect of changing interindividual distances.
include("rescaling.jl")


nothing

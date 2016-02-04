# module StaticAnalysis

# export Vec2, State, norm, unit, dot, cross, mean

"Vec2 is a 2D vector."
immutable Vec2
    x::Float64
    y::Float64
end

"A State contains the position and orientation of a particle."
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

import Base: norm, angle, dot
norm(u::Vec2) = hypot(u.x, u.y)
unit(u::Vec2) = u / norm(u)
angle(u::Vec2) = atan2(u.y, u.x)
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

include("run.jl")
include("alpha.jl")
include("plot.jl")
include("io.jl")
include("structure.jl")

# end

nothing

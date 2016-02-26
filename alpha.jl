# Compute α-shapes and distances from the edge of the α-shapes.

"An `AlphaShape` contains the different parts of an α-shape."
type AlphaShape
    outer::Vector{Vector{Int}}
    inner::Vector{Vector{Int}}
    degen::Vector{Vector{Int}}
    solo::Vector{Int}
end


"`mainshape` returns the largest outer shape as a vector of indices."
mainshape(shape::AlphaShape) = shape.outer[last(findmax(map(length, shape.outer)))]

"`mainshape` returns the largest outer shape as a vector of positions."
mainshape(shape::AlphaShape, p::Vector{Vec2}) = Vec2[p[i] for i in mainshape(shape)]


"alphashape computes the α-shape of a group of particles."
alphashape(α::Real, p::Vector{State}) = alphashape(α, [q.pos for q in p])

"alphashape computes the α-shape of a group of points."
function alphashape(α::Real, p::Vector{Vec2})
    n = length(p)
    r = abs(1/α)
    if α < 0
        incircle(x, xc, vc) = norm(x - xc) ≥ r
    elseif α == 0
        incircle(x, xc, vc) = dot(x - xc, vc) ≤ 0
    else
        incircle(x, xc, vc) = norm(x - xc) ≤ r
    end

    # find list of edges that are part of α-shape
    edges = Tuple{Int,Int}[]
    degen = Vector{Int}[]
    dist = zeros(n, n)
    @inbounds for i=1:n, j=i+1:n
        u = p[j] - p[i]
        d = norm(u)
        dist[i,j] = d

        # skip points that are too far away (or indistinct)
        0 < d ≤ 2r || continue

        # find center of circles
        h = sqrt(r^2 - (d/2)^2)
        vh = p[i] + u/2
        v = Vec2(-u.y, u.x) / d
        v0 = vh + h * v
        v1 = vh - h * v

        # skip edge if circle does not contain all points
        c0, c1 = true, true
        for k=1:n
            i != k != j || continue
            c0 &= incircle(p[k], v0, v)
            c1 &= incircle(p[k], v1, v)
            c1 || c0 || break
        end
        c1 && c0 && push!(degen, [i, j])
        c1 $ c0 && push!(edges, (i, j))
    end

    # find solitary particles
    solo = Int[]
    @inbounds for i=1:n
        if all(dist[i,i+1:end] .≥ 2r) && all(dist[1:i-1,i] .≥ 2r)
            push!(solo, i)
        end
    end

    # build list of points from list of edges
    outer = Vector{Int}[]
    @inbounds while !isempty(edges)
        e = pop!(edges)
        out = Int[e[1]]
        ix, vx = 1, e[2]
        while vx != out[1]
            push!(out, vx)
            i = findfirst(e -> vx in e, edges)
            e = edges[i]
            vx = e[1] == vx ? e[2] : e[1]
            deleteat!(edges, i)
        end
        push!(outer, out)
    end

    # distinguish inner and outer α-shapes
    m = length(outer)
    m > 0 || return AlphaShape(outer, outer, degen, solo)
    q = falses(m, m)
    @inbounds for i=1:m, j=1:m
        i != j || continue
        q[i,j] = inpolygon(p[outer[i][1]], p[outer[j]])
    end
    out = map(iseven, squeeze(sum(q, 2), 2))
    inner = outer[!out]
    outer = outer[out]
    AlphaShape(outer, inner, degen, solo)
end

"inpolygon tests wether a point is in a polygon."
function inpolygon(p::Vec2, poly::Vector{Vec2})
    wn = 0
    xd(e, p) = (e[2].x - e[1].x) * (p.y - e[1].y) - (p.x - e[1].x) * (e[2].y - e[1].y)
    push!(poly, poly[1]) # close polygon
    n = length(poly)
    for i in 1:n-1
        if poly[i].y <= p.y
            if poly[i+1].y > p.y && xd(poly[i:i+1], p) > 0
                wn += 1
            end
        else
            if poly[i+1].y <= p.y && xd(poly[i:i+1], p) < 0
                wn -= 1
            end
        end
    end
    wn != 0
end


# Distance from edge
# ------------------

"dist_from_edge computes the distance of each particle to the edge of the α-shape."
function dist_from_edge(α::Real, p::Vector{State})
    de = similar(p, Float64)
    shape = alphashape(α, p)
    @inbounds for i in eachindex(p)
        if i in shape.solo
            de[i] = NaN
            continue
        end
        dmin = Inf
        m = p[i].pos
        for kind in (shape.outer, shape.degen), v in kind, j=2:length(v)
            a, b = p[v[j-1]].pos, p[v[j]].pos
            u = b - a
            t = dot(m - a, u) / norm(u)
            d = norm(a + clamp(t, 0, 1) * u - m)
            dmin = min(dmin, d)
        end
        de[i] = dmin
    end
    de
end

function dist_from_edge(α::Real, p::Vector{State}, grid::Matrix{Vec2})
    d = Vector{Float64}(length(grid))
    shape = alphashape(α, p)
    @inbounds for i in eachindex(grid)
        dmin = Inf
        m = grid[i]
        for kind in (shape.outer, shape.degen), v in kind, j=2:length(v)
            a, b = p[v[j-1]].pos, p[v[j]].pos
            u = b - a
            t = dot(m - a, u) / norm(u)
            d0 = norm(a + clamp(t, 0, 1) * u - m)
            dmin = min(dmin, d0)
        end
        d[i] = dmin
        if any([inpolygon(m, Vec2[x.pos for x in p[v]]) for v in outer])
            d[i] = -dmin
        end
    end
    return d
end

function dist_from_edge(file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(file)

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates

    de = zeros(N, K)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]
        de[nz,k] = dist_from_edge(α, p[nz,k])
    end
    println("\r100%")
    de
end

function dist_from_edge_grid(file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(file)

    K = size(p, 2) # replicates
    nx, ny = size(detections, 1, 2)

    grid = [Vec2(x,y) for x in linspace(-50, 50, nx), y in linspace(-50, 50, ny)]

    de = zeros(nx * ny, K)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]
        de[:,k] = dist_from_edge(α, p[nz,k], grid)
    end
    println("\r100%")
    de
end


# Plotting
# --------

"`plotalpha` plots the precomputed α-shape of a group of particles."
function plotalpha(shape::AlphaShape, p::Vector{State})
    layers = Layer[]
    for u in shape.outer
        v = [u; u[1]]
        x = [q.pos.x for q in p[v]]
        y = [q.pos.y for q in p[v]]
        append!(layers, layer(x=x, y=y, Geom.point, Geom.path, Theme(default_color=colorant"red")))
    end
    for u in shape.inner
        v = [u; u[1]]
        x = [q.pos.x for q in p[v]]
        y = [q.pos.y for q in p[v]]
        append!(layers, layer(x=x, y=y, Geom.point, Geom.path, Theme(default_color=colorant"orange")))
    end
    for u in shape.degen
        v = [u; u[1]]
        x = [q.pos.x for q in p[v]]
        y = [q.pos.y for q in p[v]]
        append!(layers, layer(x=x, y=y, Geom.point, Geom.path, Theme(default_color=colorant"gray")))
    end
    x = [q.pos.x for q in p[shape.solo]]
    y = [q.pos.y for q in p[shape.solo]]
    append!(layers, layer(x=x, y=y, Geom.point, Theme(default_color=colorant"pink")))
    x = [q.pos.x for q in p]
    y = [q.pos.y for q in p]
    append!(layers, layer(x=x, y=y, Geom.point))
    plot(layers..., Coord.cartesian(fixed=true))
end


"`makevidalpha` saves a time sequence of α-shape plots."
function makevidalpha(α::Real, p::Matrix{State}, mask::BitMatrix)
    K = size(p, 2)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = find(mask[:,k])
        shape = alphashape(α, p[nz,k])
        h = plotalpha(shape, p[nz,k])
        draw(PNG(@sprintf("vid/%05d.png", k), 6inch, 6inch), h)
    end
    println("\r100%")
end

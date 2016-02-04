# α-shape related functions

"alphashape computes the α-shape of a group of particles."
function alphashape(α::Real, p::Vector{State})
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
        u = p[j].pos - p[i].pos
        d = norm(u)
        dist[i,j] = d

        # skip points that are too far away (or indistinct)
        0 < d ≤ 2r || continue

        # find center of circles
        h = sqrt(r^2 - (d/2)^2)
        vh = p[i].pos + u/2
        v = Vec2(-u.y, u.x) / d
        v0 = vh + h * v
        v1 = vh - h * v

        # skip edge if circle does not contain all points
        c0, c1 = true, true
        for k=1:n
            i != k != j || continue
            c0 &= incircle(p[k].pos, v0, v)
            c1 &= incircle(p[k].pos, v1, v)
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
    m > 0 || return outer, outer, degen, solo
    q = falses(m, m)
    @inbounds for i=1:m, j=1:m
        i != j || continue
        q[i,j] = inpolygon(p[outer[i][1]].pos, [v.pos for v in p[outer[j]]])
    end
    out = map(iseven, squeeze(sum(q, 2), 2))
    inner = outer[!out]
    outer = outer[out]
    outer, inner, degen, solo
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

"makevidalpha save a time sequence of α-shape plots"
function makevidalpha(α::Real, p::Matrix{State})
    K = size(x,2)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        n = countnz([v.pos.x for v in p[:,k]])
        outer, inner, degen, solo = alphashape(α, p[1:n,k])
        h = plotalpha(outer, inner, degen, solo, p[1:n,k])
        draw(PNG(@sprintf("vid/%05d.png", k), 6inch, 6inch), h)
    end
    println("\r100%")
end

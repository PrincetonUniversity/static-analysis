# Structure

using LightGraphs

function graph(adj::BitMatrix)
    n = size(adj, 1)
    g = DiGraph(n)
    for i in 1:n, j in 1:n
        if adj[i,j] > 0
            add_edge!(g, i, j)
        end
    end
    g
end

function run_structure(file::AbstractString)
    p, mask = h5read_particles(file)
    si = h5read(file, "social")

    # For each polarized swarm
    #   While there are links
    #     Remove weakest link
    #     For each remaining subgraph
    #       Compute principal components
    #       Save component size, weakest link value, aspect ratio, angle between 1st PC and swarm polarization in DataFrame

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates

    println("Building DataFrame…")
    cols = [
        (:Swarm,                Int),
        (:Polarization,         Float64),
        (:WeakestLink,          Float64),
        (:ComponentSize,        Int),
        (:ComponentAspectRatio, Float64),
        (:ComponentAngle,       Float64),
        (:NeighborAngle,        Float64),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = find(mask[:,k])

        # skip non-polarized swarm
        op, or = order_parameters(p[nz,k])
        op > 0.65 && or < 0.35 || continue

        net = si[nz,nz,k]
        m = mean(p[nz,k])
        for min in logspace(-2.5, -0.5, 25) #0:0.1:1
            # remove weakest links
            g = graph(net .>= min)

            for scc in strongly_connected_components(g)
                # need 3 points to define aspect ratio correctly
                length(scc) > 2 || continue

                # fit ellipse
                mg = mean(p[nz[scc],k])
                q = Float64[(j == 1 ? x.pos.x - mg.pos.x : x.pos.y - mg.pos.y) for x in p[nz[scc],k], j in 1:2]
                _, S, V = svd(q)
                θ = acos(dot(unit(m.vel), Vec2(V[1,1], V[2,1])))
                row = DataFrame(
                    Swarm                = k,
                    Polarization         = op,
                    WeakestLink          = min,
                    ComponentSize        = length(scc),
                    ComponentAspectRatio = sqrt(S[2]/S[1]),
                    ComponentAngle       = θ < pi/2 ? θ : θ - pi,
                    NeighborAngle        = 0,
                )

                # compute relative angle of neighbors
                for i in scc, j in neighbors(g, i)
                    ϕ = angle(p[nz[j],k].pos - p[nz[i],k].pos) - angle(p[nz[i],k].vel)
                    row[:NeighborAngle] = mod(ϕ + 3pi, 2pi) - pi
                    append!(df, row)
                end
            end
        end
    end
    println("\r100%")

    return df
end



using Gadfly

function plot_structure(df::AbstractDataFrame)
    idx = df[:ComponentAspectRatio] .< 1/√2

    dfab = by(df[idx,:], :WeakestLink) do d
        e = linspace(0, 160, 16*4)
        _, bin = hist(d[:ComponentSize], e)
        DataFrame(ComponentSize=(e[2:end]+e[1:end-1])/2, Density=bin ./ maximum(bin))
    end
    dfam = by(df[idx,:], :WeakestLink) do d
        DataFrame(ComponentSize=median(d[:ComponentSize]))
    end
    a = plot(
        layer(dfab, x=:WeakestLink, y=:ComponentSize, color=:Density, Geom.rectbin),
        layer(dfam, x=:WeakestLink, y=:ComponentSize, Geom.line, order=1),
        Scale.x_log10,
        Coord.cartesian(xmin=-2.5, xmax=-0.5, ymin=0, ymax=160))

    dfbb = by(df[idx,:], :WeakestLink) do d
        d[:ComponentAngle] = rad2deg(d[:ComponentAngle])
        e = linspace(-90, 90, 18*4)
        _, bin = hist(d[:ComponentAngle], e)
        DataFrame(ComponentAngle=(e[2:end]+e[1:end-1])/2, Density=bin ./ maximum(bin))
    end
    dfbm = by(df[idx,:], :WeakestLink) do d
        d[:ComponentAngle] = rad2deg(d[:ComponentAngle])
        DataFrame(
            ComponentAnglePos=median(d[d[:ComponentAngle].>0,:ComponentAngle]),
            ComponentAngleNeg=median(d[d[:ComponentAngle].<0,:ComponentAngle]))
    end
    b = plot(
        layer(dfbb, x=:WeakestLink, y=:ComponentAngle, color=:Density, Geom.rectbin),
        layer(dfbm, x=:WeakestLink, y=:ComponentAnglePos, Geom.line, order=1),
        layer(dfbm, x=:WeakestLink, y=:ComponentAngleNeg, Geom.line, order=1),
        Scale.x_log10,
        Scale.y_continuous(labels=x->@sprintf("%dº", x)),
        Coord.cartesian(xmin=-2.5, xmax=-0.5, ymin=-90, ymax=90))

    dfcb = by(df[idx,:], :WeakestLink) do d
        d[:NeighborAngle] = rad2deg(d[:NeighborAngle])
        e = linspace(-180, 180, 36*4)
        _, bin = hist(d[:NeighborAngle], e)
        DataFrame(NeighborAngle=(e[2:end]+e[1:end-1])/2, Density=bin ./ maximum(bin))
    end
    c = plot(
        layer(dfcb, x=:WeakestLink, y=:NeighborAngle, color=:Density, Geom.rectbin),
        Scale.x_log10,
        Scale.y_continuous(labels=x->@sprintf("%dº", x)),
        Guide.yticks(ticks=collect(-180:30:180)),
        Coord.cartesian(xmin=-2.5, xmax=-0.5, ymin=-180, ymax=180))

    # return a, b, c

    draw(PDF("component_size.pdf", 6inch, 4inch), a)
    draw(PDF("component_angle.pdf", 6inch, 4inch), b)
    draw(PDF("neighbor_angle.pdf", 6inch, 4inch), c)

    nothing
end

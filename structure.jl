# Structure

@step [get_pos, get_social_info, run_dist_from_edge] function run_structure(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]
    si = p.step[get_social_info]

    N = size(pos, 1) # max swarm size
    K = size(pos, 2) # replicates

    # distance from edge
    dem, msk = p.step[run_dist_from_edge]

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
        nz = find(mask[:,k] & msk[:,k])

        # skip non-polarized swarm
        op, or = order_parameters(pos[nz,k])
        op > 0.65 && or < 0.35 || continue

        net = si[nz,nz,k]
        m = mean(pos[nz,k])
        for min in logspace(-2.5, -0.5, 25) #0:0.1:1
            # remove weakest links
            g = graph(net .>= min)

            for scc in strongly_connected_components(g)
                # need 3 points to define aspect ratio correctly
                if length(scc) > 2
                    # fit ellipse
                    mg = mean(pos[nz[scc],k])
                    q = Float64[(j == 1 ? x.pos.x - mg.pos.x : x.pos.y - mg.pos.y) for x in pos[nz[scc],k], j in 1:2]
                    _, S, V = svd(q)
                    θ = acos(dot(unit(m.vel), Vec2(V[1,1], V[2,1])))
                    θ > pi/2 && (θ -= pi) # FIXME
                    ar = sqrt(S[2]/S[1])
                else
                    θ = NaN
                    ar = NaN
                end
                row = DataFrame(
                    Swarm                = k,
                    Polarization         = op,
                    WeakestLink          = min,
                    ComponentSize        = length(scc),
                    ComponentAspectRatio = ar,
                    ComponentAngle       = θ,
                    NeighborAngle        = NaN,
                )
                append!(df, row)

                # compute relative angle of neighbors
                for i in scc, j in neighbors(g, i)
                    ϕ = angle(pos[nz[j],k].pos - pos[nz[i],k].pos) - angle(pos[nz[i],k].vel)
                    row[:NeighborAngle] = mod2pi(ϕ + 3pi) - pi
                    append!(df, row)
                end
            end
        end
    end
    println("\r100%")

    df
end



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


function plot_structure(p::Vector{Project})
    for i in eachindex(p)
        plot_structure(p[i])
    end
end

@step [run_structure] function plot_structure(p::Project)
    print("  0%")

    df = copy(p.step[run_structure])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Convert angles from rad to deg
    df[:ComponentAngle] *= 180/pi
    df[:NeighborAngle] *= 180/pi

    # Filter neighbor angles
    nba = !isnan(df[:NeighborAngle])

    # Filter components with a well-defined orientation
    wdo = df[:ComponentAspectRatio] .< 1/√2

    function plot_col(df, col, edges, tag, colors, ticks, scale, ylabel)
        dfb = by(df, :WeakestLink) do d
            _, bin = hist(d[col], edges)
            DataFrame(Col=(edges[2:end]+edges[1:end-1])/2, Density=bin ./ sum(bin))
        end
        size(dfb, 1) > 0 || return
        dfm = by(df, :WeakestLink) do d
            if col == :ComponentAngle
                DataFrame(
                    MedianPos=median(d[d[col].>0,col]),
                    MedianNeg=median(d[d[col].<0,col]))
            else
                DataFrame(Median=median(d[col]))
            end
        end
        layers = layer(dfb, x=:WeakestLink, y=:Col, color=:Density, Geom.rectbin)
        if col == :ComponentSize
            append!(layers, layer(dfm, x=:WeakestLink, y=:Median, Geom.line, order=1, Theme(default_color=colorant"black")))
        elseif col == :ComponentAngle
            append!(layers, [
                layer(dfm, x=:WeakestLink, y=:MedianPos, Geom.line, order=1, Theme(default_color=colorant"black"));
                layer(dfm, x=:WeakestLink, y=:MedianNeg, Geom.line, order=1, Theme(default_color=colorant"black"))
            ])
        end
        h = plot(
            layers...,
            Coord.cartesian(xmin=-2.5, xmax=-0.5, ymin=edges[1], ymax=edges[end]),
            scale,
            Scale.x_log10,
            Scale.color_continuous(colormap=cmap, minvalue=colors[1], maxvalue=colors[2]),
            ticks,
            Guide.xlabel("Weakest link threshold"),
            Guide.ylabel(ylabel),
            Theme(; base_theme...))

        name = string(tag, "_", p.uid, ".pdf")
        draw(PDF(joinpath(plot_path, name), 6inch, 4inch), h)
    end

    print("\r 10%")

    plot_col(df[!nba,:], :ComponentSize, 0:5:160, "component_size",
        (0.0, 1.0),
        Scale.y_continuous,
        Guide.yticks,
        "Component size")
    print("\r 30%")

    plot_col(df[!nba & wdo,:], :ComponentAngle, -90:5:90, "component_angle",
        (0.0, 0.25),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)),
        Guide.yticks,
        "Component relative orientation")
    print("\r 50%")

    plot_col(df[nba,:], :NeighborAngle, -180:10:180, "neighbor_angle",
        (0.0, 0.05),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)),
        Guide.yticks(ticks=collect(-180:30:180)),
        "Relative angle of neighbors")
    print("\r100%")
end

function plot_structure_merged(p::Vector{Project})
    @assert length(p) > 0
    cmap = eval(p[1].conf[:colormap])
    base_theme = p[1].conf[:theme]
    dfm = []
    translate = Dict{AbstractString,AbstractString}(
        "random_pos_density" => "PD",
        "random_pos" => "P",
        "control" => "C",
    )
    for i in eachindex(p)
        run_structure(p[i]) # in case it hadn't run
        df = copy(p[i].step[run_structure])
        df[:Dataset] = translate[p[i].uid]
        if i == 1
            dfm = df
        else
            append!(dfm, df)
        end
    end

    print(" 30%")

    # Convert angles from rad to deg
    dfm[:ComponentAngle] *= 180/pi
    dfm[:NeighborAngle] *= 180/pi

    # Filter neighbor angles
    nba = !isnan(dfm[:NeighborAngle])

    # Filter components with a well-defined orientation
    wdo = dfm[:ComponentAspectRatio] .< 1/√2

    function plot_col(df, col, edges, tag, colors, ticks, scale, ylabel)
        dfq = by(df, [:WeakestLink, :Dataset]) do d
            if col == :ComponentAngle
                DataFrame(
                    MedianPos=median(d[d[col].>0,col]),
                    MedianNeg=median(d[d[col].<0,col]))
            else
                DataFrame(Median=median(d[col]))
            end
        end
        size(dfq, 1) > 0 || return
        layers = []
        if col == :ComponentSize
            layers = layer(dfq, x=:WeakestLink, y=:Median, color=:Dataset, Geom.line)
        elseif col == :ComponentAngle
            layers = layer(dfq, x=:WeakestLink, y=:MedianPos, color=:Dataset, Geom.line)
            append!(layers, layer(dfq, x=:WeakestLink, y=:MedianNeg, color=:Dataset, Geom.line))
        end
        h = plot(
            layers...,
            Coord.cartesian(xmin=-2.5, xmax=-0.5, ymin=edges[1], ymax=edges[end]),
            scale,
            Scale.x_log10,
            # Scale.color_continuous(colormap=cmap, minvalue=colors[1], maxvalue=colors[2]),
            ticks,
            Guide.xlabel("Weakest link threshold"),
            Guide.ylabel(ylabel),
            Theme(; base_theme...))

        name = string(tag, "_merged.pdf")
        try
            draw(PDF(joinpath(plot_path, name), 6inch, 4inch), h)
        catch
            return dfq
        end
    end

    print("\r 40%")

    plot_col(dfm[!nba,:], :ComponentSize, 0:5:160, "component_size",
        (0.0, 1.0),
        Scale.y_continuous,
        Guide.yticks,
        "Component size")
    print("\r 30%")

    plot_col(dfm[!nba & wdo,:], :ComponentAngle, -90:5:90, "component_angle",
        (0.0, 0.25),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)),
        Guide.yticks,
        "Component relative orientation")
    print("\r 50%")

    plot_col(dfm[nba,:], :NeighborAngle, -180:10:180, "neighbor_angle",
        (0.0, 0.05),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)),
        Guide.yticks(ticks=collect(-180:30:180)),
        "Relative angle of neighbors")
    print("\r100%")
end



"""
`path_weights` computes the weights of all the edges of a graph
and the multiplicative weigths of all paths of length 2.
"""
function path_weights(file::AbstractString)
    p, mask = h5read_particles(file)
    si = h5read(file, "social")

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates

    edges = 0:0.025:1
    bins = (edges[1:end-1] + edges[2:end]) / 2
    n = length(bins)

    println("Building DataFrame…")
    df = DataFrame(
        PathLength=repeat([1,2], inner=[n]),
        Weight=repeat(collect(bins), outer=[2]),
        BinCount=zeros(Int, 2n))
    @inbounds for k=1:K
        # @printf("\r%3d%%", 100(k-1) / K)
        @printf("\r%d", k)
        nz = find(mask[:,k])
        lz = length(nz)

        # skip non-polarized swarm
        op, or = order_parameters(p[nz,k])
        op > 0.65 && or < 0.35 || continue

        net = si[nz,nz,k]

        for j in 1:lz, i in 1:lz
            i != j || continue

            for l in 1:n
                if net[i,j] < edges[l+1]
                    df[l, :BinCount] += 1
                    break
                end
            end

            for m in 1:lz
                j != m || continue
                for l in 1:n
                    if net[i,j] * net[j,m] < edges[l+1]
                        df[n+l, :BinCount] += 1
                        break
                    end
                end
            end

        end
    end
    println("\r100%")

    return df
end

function plot_path_weights(df::AbstractDataFrame)
    p = plot(df, x=:Weight, y=:BinCount, color=:PathLength,
        Geom.bar(position=:dodge),
        Scale.y_log10, Scale.color_discrete)
    draw(PDF("path_weights.pdf", 6inch, 4inch), p)
end
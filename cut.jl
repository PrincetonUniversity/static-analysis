# Test cut transversal connections

function cut_transversal_connections!(net::Matrix{Float64}, pos::Vector{State})
    # for each two points, compute angle relative to dir
    # and set weight to zero if abs(angle ± 90°) < 30°
    dir = mean(pos).vel
    I, J, V = findnz(net)
    for k in eachindex(I, J)
        i, j = I[k], J[k]
        u = pos[j].pos - pos[i].pos
        θ = min(abs(angle(u, rotate(dir, π/2))), abs(angle(u, rotate(dir, -π/2))))
        if θ < π/6
            net[i,j] /= 2
        end
    end
end


@step [get_pos, get_social_info, get_detect_file, get_grid, run_dist_from_edge_grid] function run_spread_cut(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]
    si = p.step[get_social_info]
    detections = load_detections(p.step[get_detect_file])
    grid = p.step[get_grid]
    sensitivity = p.conf[:sensitivity]

    N = size(pos, 1) # max swarm size
    K = size(pos, 2) # replicates
    nx, ny = size(grid)

    gxmin, gxmax = extrema([v.x for v in grid])
    gymin, gymax = extrema([v.y for v in grid])

    # distance from edge
    dem, msk = p.step[run_dist_from_edge_grid]
    de_max = maximum(dem)

    # distances, order parameters, dataframe
    cols = [
        (:DistFromEdge,         Float64),
        (:MaxRadius,            Float64),
        (:RelativeDir,          Float64),
        (:SwarmID,              Int),
        (:Polarization,         Float64),
        (:Rotation,             Float64),
        (:State,                ASCIIString),
        (:Detections,           Int),
        (:ExpectedCascadeSize,  Float64),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]
        mz = msk[:,k]

        v = mainshape(alphashape(α, pos[nz,k])) # FIXME: slow and stupid
        xmin, xmax = extrema([pos[nz[i],k].pos.x for i in v])
        ymin, ymax = extrema([pos[nz[i],k].pos.y for i in v])
        rmax = min(abs(gxmin - xmin), abs(gxmax - xmax),
                   abs(gymin + ymin), abs(gymax - ymax))

        # order parameters and state, relative direction
        dir = NaN
        op, or = order_parameters(pos[nz,k])
        if op > 0.65 && or < 0.35
            state = "Polarized"
            dir = relative_dir(pos[nz,k], grid)[mz]
        elseif op < 0.35 && or > 0.65
            state = "Milling"
        elseif op < 0.35 && or < 0.35
            state = "Swarming"
        else
            state = "Transitioning"
        end

        net = si[nz,nz,k]

        # CUT!
        cut_transversal_connections!(net, pos[nz,k])

        # compute cascade size
        ex = expected_cascade_size(detections[:,:,k], net, nz, sensitivity)

        # append to dataframe
        row = DataFrame(
            DistFromEdge        = dem[mz,k],
            MaxRadius           = rmax,
            RelativeDir         = dir,
            SwarmID             = k,
            Polarization        = op,
            Rotation            = or,
            State               = state,
            Detections          = map(countnz, detections[:,:,k][mz]),
            ExpectedCascadeSize = ex[mz],
        )
        append!(df, row)
    end
    println("\r100%")

    df
end


@step [run_spread_cut] function plot_rel_spread_angle_dist_cut(p::Project)
    print("  0%")

    df = copy(p.step[run_spread_cut])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
    if !any(df[:State].=="Polarized")
       return
    end
    df = df[df[:State].=="Polarized",:]

    # Limit to common disk
    rmax = minimum(df[:MaxRadius])
    df = df[df[:DistFromEdge] .<= rmax,:]

    # Covert relative dir to degrees
    df[:RelativeDir] = rad2deg(mod(df[:RelativeDir] + 3pi, 2pi) - pi)

    df[:Detections] = convert(DataVector{Float64}, df[:Detections])

    #dp[dp[:Detections].==0,:Detections] = NA
    ey = -180:10:180
    for i in 1:size(df, 1)
        for j in 2:length(ey)
            if df[i,:RelativeDir] < ey[j]
                df[i,:RelativeDir] = (ey[j-1] + ey[j]) / 2
                break
            end
        end
    end

    dfd = by(df, :RelativeDir) do d
        e = 0:0.05:1
        _, bin = hist((d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections], e)
        DataFrame(RelativeSpread=(e[2:end]+e[1:end-1])/2, Density=bin ./ sum(bin))
    end
    dfm = by(df, :RelativeDir) do d
        e = 0:0.05:1
        m = median((d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections])
        DataFrame(RelativeSpread=(e[2:end]+e[1:end-1])/2, Median=m)
    end
    print("\r 50%")

    h = plot(layer(dfd, x=:RelativeDir, y=:RelativeSpread, color=:Density, Geom.rectbin),
        layer(dfm, x=:RelativeDir, y=:Median, Geom.line,
            order=1, Theme(default_color=colorant"black")),
        Coord.cartesian(xmin=-180, xmax=180),
        Scale.color_continuous(colormap=cmap),
        Scale.x_continuous(labels=x->@sprintf("%d°", x)),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Angle relative to group direction"),
        Guide.ylabel("Expected relative spread"),
        Theme(; base_theme...))

    draw(PDF(joinpath(plot_path, "rel_spread_angle_dist_cut_$(p.uid).pdf"), 6inch, 4inch), h)

    println("\r100%")
end

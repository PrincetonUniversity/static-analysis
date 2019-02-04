# Strip detections based on visual zones (binocular, blind…)

@step [get_pos, get_grid, get_detect_file, run_dist_from_edge_grid] function run_detect_blind(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]
    grid = p.step[get_grid]
    detections = load_detections(p.step[get_detect_file])
    strip_blind_detections!(detections, pos, mask, grid)

    N = size(pos, 1) # max swarm size
    K = size(pos, 2) # replicates
    nx, ny = size(grid)

    gxmin, gxmax = extrema([v.x for v in grid])
    gymin, gymax = extrema([v.y for v in grid])

    # distance from edge
    dem, msk = p.step[run_dist_from_edge_grid]

    # distances, order parameters, dataframe
    cols = [
        (:DistFromEdge,       Float64),
        (:MaxRadius,          Float64),
        (:RelativeDir,        Float64),
        (:SwarmID,            Int),
        (:Polarization,       Float64),
        (:Rotation,           Float64),
        (:State,              ASCIIString),
        (:Detections,         Int),
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

        if state != "Polarized"
            dir = relative_dir([State(Vec2(0,0), Vec2(1,0))], grid)[mz]
        end

        # append to dataframe
        row = DataFrame(
            DistFromEdge       = dem[mz,k],
            MaxRadius          = rmax,
            RelativeDir        = dir,
            SwarmID            = k,
            Polarization       = op,
            Rotation           = or,
            State              = state,
            Detections         = Vector{Int}(map(countnz, detections[:,:,k][mz])),
        )
        append!(df, row)
    end
    pool!(df, :State)
    println("\r100%")

    df
end


@step [run_detect_blind] function plot_detect_blind(p::Project)
    print("  0%")

    df = copy(p.step[run_detect_blind])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
    if all(df[:State].!="Polarized")
        return
    end
    df = df[df[:State].=="Polarized",:]

    # Limit to common disk
    rmax = minimum(df[:MaxRadius])
    df = df[df[:DistFromEdge] .<= rmax,:]

    # Covert relative dir to degrees
    df[:RelativeDir] = rad2deg(mod(df[:RelativeDir] + 3pi, 2pi) - pi)

    function bin_col(col, edges)
        for row in eachrow(df)
            for j in 2:length(edges)
                if row[col] < edges[j]
                    row[col] = (edges[j-1] + edges[j]) / 2
                    break
                end
            end
        end
        df[df[col] .< last(edges),:]
    end

    df = bin_col(:RelativeDir, -180:5:180)
    print("\r 50%")

    edges = 0:2:80
    dfd = by(df, :RelativeDir) do d
        _, bin = hist(d[:Detections], edges)
        DataFrame(Detections=(edges[2:end]+edges[1:end-1])/2, Density=bin ./ sum(bin))
    end

    h = plot(dfd, x=:RelativeDir, y=:Detections, color=:Density, Geom.rectbin,
        Coord.cartesian(xmin=-180, xmax=180),
        Scale.x_continuous(labels=x->@sprintf("%d°", x)),
        Scale.color_continuous(colormap=cmap, minvalue=0.0, maxvalue=0.15),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Relative angular position"),
        Guide.ylabel("Detection count"),
        Theme(; base_theme...))

    draw(PDF(joinpath(plot_path, "detections_angle_$(p.uid).pdf"), 6inch, 4inch), h)
    println("\r100%")
end


@step [run_detect_blind] function plot_detect_polar_blind(p::Project)
    print("  0%")

    df = copy(p.step[run_detect_blind])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
    if all(df[:State].!="Polarized")
        return
    end
    df = df[df[:State].=="Polarized",:]

    # Limit to common disk
    rmax = minimum(df[:MaxRadius])
    df = df[df[:DistFromEdge] .<= rmax,:]

    # Covert relative dir to degrees
    df[:RelativeDir] = rad2deg(mod(df[:RelativeDir] + 3pi, 2pi) - pi)

    df[:Detections] = convert(DataVector{Float64}, df[:Detections])
    ex = 0:1:30
    ey = -180:5:180
    for i in 1:size(df, 1)
        for j in 2:length(ex)
            if df[i,:DistFromEdge] < ex[j]
                df[i,:DistFromEdge] = (ex[j-1] + ex[j]) / 2
                break
            end
        end
        for j in 2:length(ey)
            if df[i,:RelativeDir] < ey[j]
                df[i,:RelativeDir] = (ey[j-1] + ey[j]) / 2
                break
            end
        end
    end
    df = df[df[:DistFromEdge].<ex[end],:]
    print("\r 50%")

    dfm = by(df, [:DistFromEdge, :RelativeDir]) do d
        DataFrame(Detections=mean(d[:Detections]))
    end

    h = plot(dfm, x=:DistFromEdge, y=:RelativeDir, color=:Detections,
        Geom.rectbin,
        Coord.cartesian(xmin=ex[1], xmax=ex[end], ymin=ey[1], ymax=ey[end]),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)),
        Scale.color_continuous(colormap=cmap, maxvalue=50),
        Guide.yticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Relative angular position"),
        Guide.colorkey("Average\ndetection\ncount"),
        Theme(; base_theme...))

    draw(PDF(joinpath(plot_path, "detections_dist_angle_$(p.uid).pdf"), 6inch, 4inch), h)

    println("\r100%")
end


@step [get_pos, get_social_info, get_detect_file, get_grid, run_dist_from_edge_grid] function run_spread_blind(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]
    si = p.step[get_social_info]
    detections = load_detections(p.step[get_detect_file])
    grid = p.step[get_grid]
    sensitivity = p.conf[:sensitivity]
    strip_blind_detections!(detections, pos, mask, grid)

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


@step [run_spread_blind] function plot_rel_spread_angle_dist_blind(p::Project)
    print("  0%")

    df = copy(p.step[run_spread_blind])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
    if all(df[:State].!="Polarized")
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
        rho = (d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections]
        _, bin = hist(rho[!isnan(rho)], e)
        DataFrame(RelativeSpread=(e[2:end]+e[1:end-1])/2, Density=bin ./ sum(bin))
    end
    dfm = by(df, :RelativeDir) do d
        e = 0:0.05:1
        rho = (d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections]
        m = median(rho[!isnan(rho)])
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

    draw(PDF(joinpath(plot_path, "rel_spread_angle_dist_blind_$(p.uid).pdf"), 6inch, 4inch), h)

    println("\r100%")
end


@step [run_spread_blind] function plot_cascade_radius_blind(p::Project)
    print("  0%")

    df = copy(p.step[run_spread_blind])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to common disk
    rmax = minimum(df[:MaxRadius])
    df = df[df[:DistFromEdge] .<= rmax,:]
    print("\r 10%")

    edges = 0:2:30
    @inbounds for row in eachrow(df)
        for j in 2:length(edges)
            if row[:DistFromEdge] < edges[j]
                row[:DistFromEdge] = (edges[j-1] + edges[j]) / 2
                break
            end
        end
    end
    df = df[df[:DistFromEdge] .< last(edges),:]
    print("\r 40%")

    dfm = by(df, [:DistFromEdge, :State]) do d
        rho = d[:ExpectedCascadeSize]
        m, s = mean(rho), sem(rho)
        DataFrame(Mean=m, Min=m-s, Max=m+s)
    end

    print("\r 60%")
    h = plot(layer(dfm, x=:DistFromEdge, ymin=:Min, ymax=:Max, color=:State, Geom.errorbar),
        layer(dfm, x=:DistFromEdge, y=:Mean, color=:State, Geom.line, order=1),
        layer(dfm, x=:DistFromEdge, y=:Mean, color=:State, Geom.point, order=2),
        Coord.cartesian(xmin=0, xmax=30, ymin=20, ymax=40),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Expected cascade size"),
        Theme(; base_theme...))

    name = joinpath(plot_path, "cascade_radius_blind_$(p.uid).pdf")
    draw(PDF(name, 6inch, 4inch), h)

    println("\r100%")
end


@step [run_spread_blind] function plot_relative_spread_radius_blind(p::Project)
    print("  0%")

    df = copy(p.step[run_spread_blind])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to common disk
    rmax = minimum(df[:MaxRadius])
    df = df[df[:DistFromEdge] .<= rmax,:]
    print("\r 10%")

    edges = 0:2:30
    @inbounds for row in eachrow(df)
        for j in 2:length(edges)
            if row[:DistFromEdge] < edges[j]
                row[:DistFromEdge] = (edges[j-1] + edges[j]) / 2
                break
            end
        end
    end
    df = df[df[:DistFromEdge] .< last(edges),:]
    print("\r 40%")

    dfm = by(df, [:DistFromEdge, :State]) do d
        rho = (d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections]
        m, s = mean(rho), sem(rho)
        DataFrame(Mean=m, Min=m-s, Max=m+s)
    end

    print("\r 60%")
    h = plot(layer(dfm, x=:DistFromEdge, ymin=:Min, ymax=:Max, color=:State, Geom.errorbar),
        layer(dfm, x=:DistFromEdge, y=:Mean, color=:State, Geom.line, order=1),
        layer(dfm, x=:DistFromEdge, y=:Mean, color=:State, Geom.point, order=2),
        Coord.cartesian(xmin=0, xmax=30, ymin=0.25, ymax=0.5),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Relative spread of cascades"),
        Theme(; base_theme...))

    name = joinpath(plot_path, "relative_spread_radius_blind_$(p.uid).pdf")
    draw(PDF(name, 6inch, 4inch), h)

    println("\r100%")
end


function plot_avg_relative_spread_radius_blind(p::Vector{Project})
    print("  0%")

    cmap = eval(p[1].conf[:colormap])
    base_theme = p[1].conf[:theme]

    dfm = DataFrame([Float64, Float64, Float64, Float64, Float64, ASCIIString], [:Q1, :Q2, :Q3, :Top, :Bottom, :Mode], 0)

    mapping = Dict{AbstractString, AbstractString}(
        "random_pos_vel_density" => "POD",
        "random_pos_vel" => "PO",
        "random_pos_density" => "PD",
        "random_pos" => "P",
        "random_vel" => "O",
        "control" => "C",
    )

    for i in eachindex(p)
        df = p[i].step[run_spread_blind]
        rho = (df[:ExpectedCascadeSize] - df[:Detections]) ./ df[:Detections]
        if !all(isnan(rho))
            rho = rho[!isnan(rho)]
            Q1, Q2, Q3 = percentile(rho, 25), median(rho), percentile(rho, 75)
            IQR = Q3 - Q1
            sorted = sort(rho)
            top = sorted[sorted.<=(Q3+1.5*IQR)][end]
            bot = sorted[sorted.>=(Q1-1.5*IQR)][1]
            # outliers = sorted[sorted.>(Q3+1.5*IQR)]
            # append!(outliers, sorted[sorted.<(Q1-1.5*IQR)])
            df = DataFrame(Q1=Q1, Q2=Q2, Q3=Q3, Top=top, Bottom=bot, Mode=mapping[p[i].uid])
            append!(dfm, df)
        end
    end

    print("\r 60%")
    h = plot(dfm, x=:Mode, middle=:Q2, lower_hinge=:Q1, upper_hinge=:Q3,
        lower_fence=:Bottom, upper_fence=:Top,
        Geom.boxplot, Stat.identity,
        Coord.cartesian(ymin=0, ymax=0.9),
        Guide.xlabel("Random mode"),
        Guide.ylabel("Avg. rel. spread of cascades"),
        Theme(; base_theme...))

    name = joinpath(plot_path, "avg_relative_spread_radius_blind.pdf")
    draw(PDF(name, 6inch, 4inch), h)

    println("\r100%")
end


function plot_avg_relative_spread_radius_blind_vs_not(p::Vector{Project})
    print("  0%")

    cmap = eval(p[1].conf[:colormap])
    base_theme = p[1].conf[:theme]

    dfm = DataFrame([Float64, Float64, Float64, Float64, Float64, ASCIIString, UTF8String], [:Q1, :Q2, :Q3, :Top, :Bottom, :Mode, :BlindAngle], 0)

    mapping = Dict{AbstractString, AbstractString}(
        "random_pos_vel_density" => "POD",
        "random_pos_vel" => "PO",
        "random_pos_density" => "PD",
        "random_pos" => "P",
        "random_vel" => "O",
        "control" => "C",
    )

    for i in eachindex(p), step in [run_spread, run_spread_blind]
        df = p[i].step[step]
        rho = (df[:ExpectedCascadeSize] - df[:Detections]) ./ df[:Detections]
        if !all(isnan(rho))
            rho = rho[!isnan(rho)]
            Q1, Q2, Q3 = percentile(rho, 25), median(rho), percentile(rho, 75)
            IQR = Q3 - Q1
            sorted = sort(rho)
            top = sorted[sorted.<=(Q3+1.5*IQR)][end]
            bot = sorted[sorted.>=(Q1-1.5*IQR)][1]
            # outliers = sorted[sorted.>(Q3+1.5*IQR)]
            # append!(outliers, sorted[sorted.<(Q1-1.5*IQR)])
            df = DataFrame(Q1=Q1, Q2=Q2, Q3=Q3, Top=top, Bottom=bot, Mode=mapping[p[i].uid], BlindAngle=(step==run_spread?"0°":"25°"))
            append!(dfm, df)
        end
    end
    dfm
end

function plot_avg_relative_spread_radius_blind_vs_not2(p::Vector{Project}, df::AbstractDataFrame)
    cmap = eval(p[1].conf[:colormap])
    base_theme = p[1].conf[:theme]

    dfm = deepcopy(df)
    for row in eachrow(dfm)
        row[:Mode] = string(row[:Mode], row[:BlindAngle]=="0°"?"1":"2")
    end

    print("\r 60%")
    h = plot(dfm, x=:Mode, middle=:Q2, lower_hinge=:Q1, upper_hinge=:Q3,
        lower_fence=:Bottom, upper_fence=:Top, color=:BlindAngle,
        Geom.boxplot, Stat.identity,
        Coord.cartesian(ymin=0, ymax=0.9),
        Guide.xlabel("Random mode"),
        Guide.ylabel("Avg. rel. spread of cascades"),
        Guide.colorkey("Blind angle"),
        Guide.xticks(ticks=[1:2:11;]),
        Scale.x_discrete(labels=x->strip(x, '1')),
        Theme(; base_theme...))

    name = joinpath(plot_path, "avg_relative_spread_radius_blind_vs_not.pdf")
    draw(PDF(name, 6inch, 4inch), h)

    println("\r100%")
end



# function load_detections(detect_file::AbstractString)
#     p, mask = h5read_particles(detect_file)
#     detections = h5read(detect_file, "detections")
#     p, mask, detections
# end

function strip_non_binocular_detections!(detections::Array{BitVector,3},
                                         p::Matrix{State}, mask::BitMatrix,
                                         grid::Matrix{Vec2},
                                         binocular_angle::Real=deg2rad(31/2))
    @assert(0 ≤ binocular_angle ≤ π)
    strip_detections!(detections, p, mask, grid, π - binocular_angle)
end

function strip_blind_detections!(detections::Array{BitVector,3},
                                 p::Matrix{State}, mask::BitMatrix,
                                 grid::Matrix{Vec2},
                                 blind_angle::Real=deg2rad(25/2))
    @assert(0 ≤ blind_angle ≤ π)
    strip_detections!(detections, p, mask, grid, blind_angle)
end

"""
`strip_detections!` removes detections located within an angle
`angle_from_back` from the back of each particle.
"""
function strip_detections!(detections::Array{BitVector,3},
                           p::Matrix{State}, mask::BitMatrix,
                           grid::Matrix{Vec2},
                           angle_from_back::Real)
    K = size(p, 2)
    nx, ny = size(detections, 1, 2)
    for k in 1:K, j in 1:nx, i in 1:ny
        for l in find(detections[i,j,k])
            if abs(angle(p[l,k].pos - grid[i,j], p[l,k].vel)) ≤ angle_from_back
                detections[i,j,k][l] = false
            end
        end
    end
end


# Plotting
# --------

function plot_zones(p::Matrix{State}, mask::BitMatrix,
                    detections_all::Array{BitVector, 3},
                    detections_binoc::Array{BitVector, 3},
                    detections_blind::Array{BitVector, 3},
                    grid::Matrix{Vec2})
    k = 1
    nz = find(mask[:,k])

    draw(PDF("detect_all.pdf", 6inch, 4inch),
        plotdetect(p[nz,k], detections_all[:,:,k], grid))
    draw(PDF("detect_binoc.pdf", 6inch, 4inch),
        plotdetect(p[nz,k], detections_binoc[:,:,k], grid))
    draw(PDF("detect_blind.pdf", 6inch, 4inch),
        plotdetect(p[nz,k], detections_blind[:,:,k], grid))
end

function plot_zones_handpicked(p::Matrix{State}, mask::BitMatrix,
                               detections_all::Array{BitVector, 3},
                               detections_binoc::Array{BitVector, 3},
                               detections_blind::Array{BitVector, 3},
                               grid::Matrix{Vec2})
    for (k, name) in enumerate(["polarized", "milling2", "milling", "swarming"])
        nz = find(mask[:,k])
        m = maximum(map(countnz, detections_all[:,:,k]))
        draw(PDF("$(name)_all.pdf", 5inch, 5inch),
            plotdetect(p[nz,k], detections_all[:,:,k], grid, minvalue=0, maxvalue=m, title="full vision — $(name)"))
        draw(PDF("$(name)_binoc.pdf", 5inch, 5inch),
            plotdetect(p[nz,k], detections_binoc[:,:,k], grid, minvalue=0, title="binocular only — $(name)"))
        draw(PDF("$(name)_blind.pdf", 5inch, 5inch),
            plotdetect(p[nz,k], detections_blind[:,:,k], grid, minvalue=0, maxvalue=m, title="with blind zone — $(name)"))
    end
end
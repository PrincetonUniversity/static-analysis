# Analyze detection count in and around the schools.

@step [get_detect_file] function get_grid(p::Project)
    load_grid(p.step[get_detect_file])
end

# @step [get_detect_file] function get_detections(p::Project)
#     load_detections(p.step[get_detect_file])
# end

@step [get_pos, get_grid] function run_dist_from_edge_grid(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]

    K = size(pos, 2) # replicates

    grid = p.step[get_grid]
    nx, ny = size(grid)

    de = zeros(nx * ny, K)
    msk = falses(nx * ny, K)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]
        de[:,k], msk[:,k] = dist_from_edge(α, pos[nz,k], grid)
    end
    println("\r100%")
    de, msk
end


@step [get_pos, get_grid, get_detect_file, run_dist_from_edge_grid] function run_detect(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]
    grid = p.step[get_grid]
    detections = load_detections(p.step[get_detect_file])

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


@step [run_detect] function plot_detect(p::Project)
    print("  0%")

    df = copy(p.step[run_detect])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
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
        Scale.color_continuous(colormap=cmap),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Relative angular position"),
        Guide.ylabel("Detection count"),
        Theme(; base_theme...))

    draw(PDF(joinpath(plot_path, "detections_angle.pdf"), 6inch, 4inch), h)
    println("\r100%")
end


@step [run_detect] function plot_detect_polar(p::Project)
    print("  0%")

    df = copy(p.step[run_detect])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
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

    draw(PDF(joinpath(plot_path, "detections_dist_angle.pdf"), 6inch, 4inch), h)

    println("\r100%")
end


@step [run_detect] function plot_detect_radius(p::Project)
    print("  0%")

    df = copy(p.step[run_detect])
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
        m, s = mean(d[:Detections]), sem(d[:Detections])
        DataFrame(Mean=m, Min=m-s, Max=m+s)
    end

    print("\r 60%")
    h = plot(layer(dfm, x=:DistFromEdge, ymin=:Min, ymax=:Max, color=:State, Geom.errorbar),
        layer(dfm, x=:DistFromEdge, y=:Mean, color=:State, Geom.line, order=1),
        layer(dfm, x=:DistFromEdge, y=:Mean, color=:State, Geom.point, order=2),
        Coord.cartesian(xmin=0, xmax=30, ymin=15, ymax=35),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Number of detections"),
        Theme(; base_theme...))

    name = joinpath(plot_path, "detect_radius.pdf")
    draw(PDF(name, 6inch, 4inch), h)

    println("\r100%")
end


# I/O
# ---

"`load_detections` loads detections from `file` as an array of `BitVector`s."
function load_detections(file::AbstractString)
    detections = h5read(file, "detections")
    convert_detections(detections)
end

function convert_detections(detections::Array{Array{UInt64},3})
    function makeBitVector(b::Vector{UInt64})
        v = falses(256)
        v.chunks = b
        v
    end
    map(makeBitVector, detections)
end

function load_grid(file::AbstractString)
    conf = h5readattr(file, "config")
    [Vec2(x,y) for
        x in linspace(conf["GridXmin"], conf["GridXmax"], conf["GridXcount"]),
        y in linspace(conf["GridYmin"], conf["GridYmax"], conf["GridYcount"])]
end


# Processing
# ----------

function run_detect(detections::Array{BitVector,3},
                    p::Matrix{State}, mask::BitMatrix,
                    grid::Matrix{Vec2};
                    α::Real = -0.2)

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates
    nx, ny = size(detections, 1, 2)

    @assert(size(detections, 3) == K)
    @assert(size(grid) == (nx, ny))

    # distance from edge
    println("Loading distance from edge…")
    dem = h5read(joinpath(data_path, "dist_from_edge_grid.h5"), "DistFromEdge")
    de_max = maximum(dem)

    # distances, order parameters, dataframe
    println("Building DataFrame…")
    cols = [
        (:DistFromEdge,       Float64),
        (:InsideMaxDisk,      Bool),
        (:RelativeDir,        Float64),
        (:Swarm,              Int),
        (:SwarmPolarization,  Float64),
        (:SwarmRotation,      Float64),
        (:SwarmState,         Symbol),
        (:Detections,         Int),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]

        xmin, xmax = extrema([s.pos.x for s in p[nz,k]])
        ymin, ymax = extrema([s.pos.y for s in p[nz,k]])
        max = min(50 + xmin, 50 - xmax, 50 + ymin, 50 - ymax)

        # order parameters and state, relative direction
        dir = NaN
        op, or = order_parameters(p[nz,k])
        if op > 0.65 && or < 0.35
            state = :Polarized
            dir = relative_dir(p[nz,k], grid)
        elseif op < 0.35 && or > 0.65
            state = :Milling
        elseif op < 0.35 && or < 0.35
            state = :Swarming
        else
            state = :Transitioning
        end

        # append to dataframe
        row = DataFrame(
            DistFromEdge       = dem[:,k],
            InsideMaxDisk      = dem[:,k] .<= max,
            RelativeDir        = dir,
            Swarm              = k,
            SwarmPolarization  = op,
            SwarmRotation      = or,
            SwarmState         = state,
            Detections         = Vector{Int}(map(countnz, detections[:,:,k][:])),
        )
        append!(df, row)
    end
    println("\r100%")

    return df
end


function relative_dir(p::Vector{State}, grid::Matrix{Vec2})
    m = mean(p)
    θ = angle(m.vel)
    [angle(grid[i] - m.pos) - θ for i in eachindex(grid)]
end


"dist_from_edge computes the distance of each grid point to the edge of the α-shape."
function dist_from_edge(α::Real, p::Vector{State}, grid::Matrix{Vec2})
    mask = falses(length(grid))
    de = Vector{Float64}(length(grid))
    shape = alphashape(α, p)
    v = mainshape(shape)
    vp = Vec2[p[i].pos for i in v]
    @inbounds for i in eachindex(grid)
        m = grid[i]
        dmin = Inf
        for j=1:length(v)
            k = j == 1 ? length(v) : j-1
            a, b = p[v[j]].pos, p[v[k]].pos
            u = b - a
            t = dot(m - a, u) / norm(u)^2
            d = norm(a + clamp(t, 0, 1) * u - m)
            dmin = min(dmin, d)
        end
        de[i] = dmin
        if inpolygon(m, vp)
            de[i] = -dmin
        else
            mask[i] = true
        end
    end
    de, mask
end


function dist_from_edge_grid(file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(file)

    K = size(p, 2) # replicates

    grid = load_grid(file)
    nx, ny = size(grid)

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

function plotgrid{T<:Real}(p::Vector{State}, val::AbstractMatrix{T}, grid::Matrix{Vec2};
                  minvalue=nothing, maxvalue=nothing, title=nothing,
                  colormap=x -> RGB(x, x, 0.5 - x/2))
    nx, ny = size(grid)
    x = [z.x for z in grid[:]]
    y = [z.y for z in grid[:]]
    df = DataFrame(X=x, Y=y, Value=val[:])
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context(units=box)
    for i in eachindex(p)
        ctx = context(rotation=Rotation(-angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
        c = compose(c, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2)))
        # c = compose(c, myellipse(p[i].pos, Vec2(1, 0.125), 0.8, angle(p[i].vel)))
    end
    plot(df, x=:X, y=:Y, color=:Value, Geom.rectbin,
        Guide.annotation(c),
        Guide.xlabel("x (body length)"),
        Guide.ylabel("y (body length)"),
        Guide.title(title),
        Scale.color_continuous(minvalue=minvalue, maxvalue=maxvalue, colormap=colormap),
        Coord.cartesian(raster=true, xmin=-50, xmax=50, ymin=-50, ymax=50, fixed=true))
end


function plotdetect(p::Vector{State}, detect::Matrix, grid::Matrix{Vec2};
                    minvalue=nothing, maxvalue=nothing, title=nothing)
    nx, ny = size(detect)
    x = [z.x for z in grid[:]]
    y = [z.y for z in grid[:]]
    df = DataFrame(X=x, Y=y, Count=map(countnz, detect[:]))
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context(units=box)
    for i in eachindex(p)
        c = compose(c, myellipse(p[i].pos, Vec2(1, 0.125), 0.8, angle(p[i].vel)))
    end
    plot(df, x=:X, y=:Y, color=:Count, Geom.rectbin,
        Guide.annotation(c),
        Guide.xlabel("x (body length)"),
        Guide.ylabel("y (body length)"),
        Guide.title(title),
        Scale.color_continuous(minvalue=minvalue, maxvalue=maxvalue, colormap=x -> RGB(x, x, 0.5 - x/2)),
        Coord.cartesian(raster=true, xmin=-50, xmax=50, ymin=-50, ymax=50, aspect_ratio=1))
end

function imagedetect(p::Vector{State}, detections::Matrix, grid::Matrix{Vec2})
    nx, ny = size(grid)
    xmin, xmax = extrema([z.x for z in grid[:]])
    ymin, ymax = extrema([z.y for z in grid[:]])
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context(units=box)
    detect = map(countnz, detections)
    m = maximum(detect)
    x, y = linspace(-50, 50, nx), linspace(-50, 50, ny)
    dx, dy = 100/nx, 100/ny
    for j in 1:ny, i in 1:nx
        c = compose(c, (context(), rectangle(grid[i,j].x-dx/2, grid[i,j].y-dy/2, 1.25dx, 1.25dy), fill(RGB(detect[i,j]/m, detect[i,j]/m, 0.5 - detect[i,j]/2m)), stroke(nothing)))
    end
    cp = context(units=box)
    for i in eachindex(p)
        cp = compose(cp, myellipse(p[i].pos, Vec2(1, 0.125), 0.8, angle(p[i].vel)))
    end
    compose(context(), cp, c)
end

function detectvid(p::Matrix{State}, mask::BitArray{2}, detect::Array{Int, 3})
    dir = "plots/vid/detect2"
    dir = joinpath(plot_path, "vid", "detect")
    try rm(dir, recursive=true) end
    try mkdir(dir) end
    K = size(p, 2)
    for k in 1:K
        @printf("\r%3d%%", 100(k-1) / K)
        sel = 1:(findfirst(!mask[:,k])-1)
        h = imagedetect(p[mask[:,k],k], detect[:,:,k])
        file = @sprintf("%06d.png", k)
        draw(PNG(joinpath(dir, file), 1000px, 1000px), h)
    end
    println("\r100%")
end

function plot_detect_radius_old(df::AbstractDataFrame)
    # Limit to max disk and outside swarm
    dp = df[df[:InsideMaxDisk] & (df[:DistFromEdge].>=0),:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])

    plot(dp, x=:DistFromEdge, y=:Detections,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Distance from boundary (body length)"),
        Guide.ylabel("Number of overlapping fields of view"),
        Guide.colorkey("State"))
end

function plot_detect_polar_alt(df::AbstractDataFrame)
    # Limit to polarized, max disk and outside swarm
    idx = (df[:SwarmState].==:Polarized) & df[:InsideMaxDisk] & (df[:DistFromEdge].>=0)

    dp = df[idx,:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])
    dp[:RelativeDir] = rad2deg(mod(dp[:RelativeDir] + 3pi, 2pi) - pi)

    dp[:BinDistFromEdge] = convert(DataVector{Float64}, zeros(size(dp, 1)))
    dp[:BinRelativeDir] = convert(DataVector{Float64}, zeros(size(dp, 1)))
    ex = 0:1:30
    ey = -180:5:180
    dp = dp[dp[:DistFromEdge].<ex[end],:]
    for i in 1:size(dp, 1)
        for j in 2:length(ex)
            if dp[i,:DistFromEdge] < ex[j]
                dp[i,:BinDistFromEdge] = (ex[j-1] + ex[j]) / 2
                break
            end
        end
        for j in 2:length(ey)
            if dp[i,:RelativeDir] < ey[j]
                dp[i,:BinRelativeDir] = (ey[j-1] + ey[j]) / 2
                break
            end
        end
    end

    dq = by(dp, [:BinDistFromEdge, :BinRelativeDir]) do d
        DataFrame(Detections=mean(d[:Detections]))
    end

    p = plot(dq, x=:BinDistFromEdge, y=:BinRelativeDir, color=:Detections,
        Geom.rectbin,
        Coord.cartesian(xmin=0, xmax=ex[end], ymin=-180, ymax=180),
        Scale.y_continuous(labels=x->@sprintf("%dº", x)),
        Scale.color_continuous(colormap=viridis, maxvalue=50),
        Guide.yticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Relative angular position"),
        Guide.colorkey("Average\ndetection\ncount"))

    draw(PDF(joinpath(plot_path, "detections_dist_angle.pdf"), 6inch, 4inch), p)
end

function plot_detect_polar2(df::AbstractDataFrame)
    # Limit to polarized, max disk and outside swarm
    idx = (df[:SwarmState].==:Polarized) & df[:InsideMaxDisk] & (df[:DistFromEdge].>=0)

    dp = df[idx,:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])
    dp[:RelativeDir] = rad2deg(mod(dp[:RelativeDir] + 3pi, 2pi) - pi)

    p = plot(dp, x=:RelativeDir, y=:Detections, Geom.line,
        Stat.binmean(n=36*2),
        Coord.cartesian(xmin=-180, xmax=180),
        Scale.x_continuous(labels=x->@sprintf("%dº", x)),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Relative angular position"),
        Guide.ylabel("Detection count"))

    draw(PDF(joinpath(plot_path, "detections_angle.pdf"), 6inch, 4inch), p)
end

function plot_detect_polar3(df::AbstractDataFrame)
    # Limit to polarized, max disk and outside swarm
    idx = (df[:SwarmState].==:Polarized) & df[:InsideMaxDisk] & (df[:DistFromEdge].>=0)

    dp = df[idx,:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])
    dp[:RelativeDir] = rad2deg(mod(dp[:RelativeDir] + 3pi, 2pi) - pi)

    dp[:BinRelativeDir] = convert(DataVector{Float64}, zeros(size(dp, 1)))
    ey = -180:5:180
    for i in 1:size(dp, 1)
        for j in 2:length(ey)
            if dp[i,:RelativeDir] < ey[j]
                dp[i,:BinRelativeDir] = (ey[j-1] + ey[j]) / 2
                break
            end
        end
    end

    dq = by(dp, :BinRelativeDir) do d
        e = 0:2:80
        _, bin = hist(d[:Detections], e)
        DataFrame(Detections=(e[2:end]+e[1:end-1])/2, Density=bin ./ maximum(bin))
    end

    p = plot(dq, x=:BinRelativeDir, y=:Detections, color=:Density, Geom.rectbin,
        Coord.cartesian(xmin=-180, xmax=180),
        Scale.x_continuous(labels=x->@sprintf("%dº", x)),
        Scale.color_continuous(colormap=viridis),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Relative angular position"),
        Guide.ylabel("Detection count"))

    draw(PDF(joinpath(plot_path, "detections_angle.pdf"), 6inch, 4inch), p)
end


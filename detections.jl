# Analyze detection count in and around the schools.


# BitSet
function unpack(b::Vector{UInt64})
    id = Int[]
    @inbounds for n in 1:256
        i = (((n-1) & 0xc0) >> 6) + 1
        j = UInt64(1) << ((n-1) & 0x3f)
        if b[i] & j > 0
            push!(id, n)
        end
    end
    id
end

bitcount(v::Vector{UInt64}) = sum(map(bitcount, v))
function bitcount(i::UInt64)
    i = i - ((i >> 1) & 0x5555555555555555);
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
    return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;
end

function set(b::Vector{UInt64}, n::Int)
    b[(n & 0xc0) >> 6] |= 1 << (n & 0x3f)
end

function unset(b::Vector{UInt64}, n::Int)
    b[(n & 0xc0) >> 6] $= ~(1 << (n & 0x3f))
end



function run_detect(file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(file)
    detections = h5read(file, "detections")

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates
    nx, ny = size(detections, 1, 2)

    # distance from edge
    grid = [Vec2(x,y) for x in linspace(-50, 50, nx), y in linspace(-50, 50, ny)]
    println("Loading distance from edge…")
    dem = h5read("$(dirname(file))/dist_from_edge_grid.h5", "DistFromEdge")
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
            Detections         = Vector{Int}(map(bitcount, detections[:,:,k][:])),
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


# Plotting
# --------

function imagedetect(p::Vector{State}, detect::Matrix)
    nx, ny = size(detect)
    box = UnitBox(-50, 50, 100, -100)
    c = context(units=box)
    m = maximum(detect)
    x, y = linspace(-50, 50, nx), linspace(-50, 50, ny)
    dx, dy = 100/nx, 100/ny
    for i in 1:nx, j in 1:ny
        c = compose(c, (context(), rectangle(x[i]-dx/2, y[j]-dy/2, 2dx, 2dy), fill(RGB(detect[i,j]/m, detect[i,j]/m, 0.5 - detect[i,j]/2m)), stroke(nothing)))
    end
    cp = context(units=box)
    for i in eachindex(p)
        cp = compose(cp, myellipse(p[i].pos, Vec2(1, 0.125), 0.8, atan2(p[i].vel.y, p[i].vel.x)))
    end
    compose(context(), cp, c)
end

function detectvid(p::Matrix{State}, mask::BitArray{2}, detect::Array{Int, 3})
    dir = "plots/vid/detect2"
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

function plotdetect(p::Vector{State}, detect::Matrix)
    nx, ny = size(detect)
    df = DataFrame(
        X=repeat(collect(linspace(-50, 50, nx)), outer=[ny]),
        Y=repeat(collect(linspace(-50, 50, ny)), inner=[nx]),
        Count=detect[:])
    box = UnitBox(-50, 50, 100, -100)
    c = context(units=box)
    for i in eachindex(p)
        c = compose(c, myellipse(p[i].pos, Vec2(1, 0.125), 0.8, atan2(p[i].vel.y, p[i].vel.x)))
    end
    plot(df, x=:X, y=:Y, color=:Count, Geom.rectbin,
        Guide.annotation(c),
        Guide.xlabel("x (body length)"),
        Guide.ylabel("y (body length)"),
        Scale.ContinuousColorScale(x -> RGB(x, x, 0.5 - x/2)),
        Coord.cartesian(raster=true, xmin=-50, xmax=50, ymin=-50, ymax=50, aspect_ratio=1))
end


function plot_detect_radius(df::AbstractDataFrame)
    # Limit to max disk and outside swarm
    dp = df[df[:InsideMaxDisk] & (df[:DistFromEdge].>=0),:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])

    plot(dp, x=:DistFromEdge, y=:Detections,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Distance from boundary (body length)"),
        Guide.ylabel("Number of overlapping fields of view"),
        Guide.colorkey("State"))
end

function plot_detect_polar(df::AbstractDataFrame)
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
        Guide.yticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Angle from school orientation"),
        Guide.title("Polarized schools only"),
        Guide.colorkey("Average\ndetection\ncount"))

    draw(PDF("detections_dist_angle.pdf", 6inch, 4inch), p)
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
        Guide.xlabel("Angle relative to group direction"),
        Guide.ylabel("Detection count"),
        Guide.title("Polarized schools only"))

    draw(PDF("detections_angle.pdf", 6inch, 4inch), p)
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
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Angle relative to group direction"),
        Guide.ylabel("Detection count"),
        Guide.title("Polarized schools only"))

    draw(PDF("detections_angle.pdf", 6inch, 4inch), p)
end


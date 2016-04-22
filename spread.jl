# Analyse correlation of information in detections.

"""
`run_spread` computes for each point in the grid the expected cascade size.
"""
function run_spread(detect_file::AbstractString, net_file::AbstractString; α::Real = -0.2)
    p, mask = h5read_particles(detect_file)
    detections = h5read(detect_file, "detections")
    si = h5read(net_file, "social")
    grid = load_grid(detect_file)

    run_spread(detections, si, p, mask, grid, α=α)
end

function run_spread(detections::Array{BitVector,3},
                    si::Array{Float64,3},
                    p::Matrix{State}, mask::BitMatrix,
                    grid::Matrix{Vec2};
                    α::Real = -0.2)
    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates
    nx, ny = size(grid)

    # distance from edge
    println("Loading distance from edge…")
    dem = h5read(joinpath(data_path, "dist_from_edge_grid_sub.h5"), "DistFromEdge")
    de_max = maximum(dem)

    # distances, order parameters, dataframe
    println("Building DataFrame…")
    cols = [
        (:DistFromEdge,        Float64),
        (:InsideMaxDisk,       Bool),
        (:RelativeDir,         Float64),
        (:Swarm,               Int),
        (:SwarmPolarization,   Float64),
        (:SwarmRotation,       Float64),
        (:SwarmState,          Symbol),
        (:Detections,          Int),
        (:ExpectedCascadeSize, Float64),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    for k=1:K
        # @printf("\r%3d%%", 100(k-1) / K)
        println(k)
        nz = find(mask[:,k])

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

        net = si[nz,nz,k]

        # compute cascade size
        ex = expected_cascade_size(detections[:,:,k], net, mask[:,k])

        # append to dataframe
        row = DataFrame(
            DistFromEdge        = dem[:,k],
            InsideMaxDisk       = dem[:,k] .<= max,
            RelativeDir         = dir,
            Swarm               = k,
            SwarmPolarization   = op,
            SwarmRotation       = or,
            SwarmState          = state,
            Detections          = map(countnz, detections[:,:,k][:]),
            ExpectedCascadeSize = ex[:],
        )
        append!(df, row)
    end
    println("\r100%")

    return df
end

function run_cascade()
    detect_file = joinpath(data_path, "handpickeddetectid.h5")
    net_file = joinpath(data_path, "handpickedprocessed.h5")
    p, mask = h5read_particles(detect_file)
    detections = load_detections(detect_file)
    grid = load_grid(detect_file)
    si = h5read(net_file, "social")

    strip_blind_detections!(detections, p, mask, grid)

    bex = similar(detections, Float64)
    label = ["polarized", "milling2", "milling", "swarming"]
    for k in 1:4
        println("# ", label[k])

        m = mask[:,k]
        ex = expected_cascade_size(detections[:,:,k], si[m,m,k], m)
        dc = map(countnz, detections[:,:,k])

        bex[:,:,k] = ex

        a = plotgrid(p[m,k], dc, grid, title="Detections")
        b = plotgrid(p[m,k], ex, grid, title="Cascade size")
        c = plotgrid(p[m,k], ex - dc, grid, title="Spread")
        d = plotgrid(p[m,k], (ex - dc) ./ dc, grid, title="Relative spread")

        path = joinpath(plot_path, "spread_blind")
        try mkdir(path) end
        draw(PDF(joinpath(path, "$(label[k])_detections.pdf"), 6inch, 6inch), a)
        draw(PDF(joinpath(path, "$(label[k])_cascade_size.pdf"), 6inch, 6inch), b)
        draw(PDF(joinpath(path, "$(label[k])_spread.pdf"), 6inch, 6inch), c)
        draw(PDF(joinpath(path, "$(label[k])_relative_spread.pdf"), 6inch, 6inch), d)
    end

    bex
end


function run_cascade_one()
    detect_file = joinpath(data_path, "handpickeddetectid.h5")
    net_file = joinpath(data_path, "handpickedprocessed.h5")
    p, mask = h5read_particles(detect_file)
    detections = load_detections(detect_file)
    grid = load_grid(detect_file)
    si = h5read(net_file, "social")

    bex = similar(detections, Float64)
    label = ["polarized", "milling2", "milling", "swarming"]
    k = 1
    println("# ", label[k])
    m = mask[:,k]
    expected_cascade_size(detections[:,:,k], si[m,m,k], m)
end


function expected_cascade_size(detect::Matrix{BitVector}, net::Matrix{Float64}, mask::BitVector)
    I, J, V = findnz(net)
    order = sortperm(V, rev=true)
    edges = zip(I[order], J[order])
    n = length(mask)

    visited = BitMatrix(size(net)...)
    r1 = BitVector(n)
    r2 = BitVector(n)

    # ex = similar(detect, Float64)
    ex = SharedArray(Float64, size(detect))
    # @inbounds for i in eachindex(ex)
    @inbounds @sync @parallel for i in eachindex(ex)
        r0 = detect[i][1:n] & mask
        ex[i] = expected_cascade_size_brute_force(r0, r1, r2, net, edges, visited)
    end
    fetch(ex)
end

function expected_cascade_size_brute_force(r0::BitVector, r1::BitVector, r2::BitVector,
                                           net::Matrix{Float64},
                                           edges, visited)
    N = 5
    n = 0.0

    @inbounds for k in 1:N
        copy!(r1, r0)
        n += expected_cascade_size_brute_force(r1, r2, net, edges, visited)
    end
    n / N
end

function expected_cascade_size_brute_force(r0::BitVector, r1::BitVector,
                                           net::Matrix{Float64},
                                           edges, visited)
    fill!(visited, false)
    expanding = true
    @inbounds while expanding
        expanding = false
        copy!(r1, r0)
        for (i, j) in edges
            if !visited[i,j] && !r0[i] && r0[j]
                visited[i,j] = true
                expanding = true
                if rand() < net[i,j]
                    r1[i] = true
                end
            end
        end
        r0 = r1
    end
    countnz(r0)
end



function expected_cascade_size_brute_force_old(detection, net, nz)
    N = 1
    n = 0.0

    I, J, V = findnz(net)
    order = sortperm(V, rev=true)
    edges = zip(I[order], J[order])
    visited = BitMatrix(size(net)...)

    @inbounds for k in 1:N
        r0 = intersect(find(detection), nz)
        visited[:] = false
        n += expected_cascade_size_brute_force_old(r0, visited, net, edges)
    end
    n / N
end

function expected_cascade_size_brute_force_old(r0::Vector{Int}, visited::BitMatrix, net::Matrix{Float64}, edges)
    oldvis = -1
    newvis = 0
    @inbounds while newvis > oldvis
        oldvis = newvis
        r1 = deepcopy(r0)
        for (i, j) in edges
            if !visited[i,j] && !(i in r0) && (j in r0)
                visited[i,j] = true
                newvis += 1
                if rand() < net[i,j]
                    push!(r1, i)
                end
            end
        end
        r0 = r1
    end
    length(r0)
end

function cart_prod(v::Vector{Float64}, d::Int)
    isempty(v) && return 0.0
    d == 1 && return sum(v)
    v[1] * cart_prod(v[2:end], d - 1) + cart_prod(v[2:end], d)
end


function plot_spread_dist(df::AbstractDataFrame; dir="spread")

    dir = joinpath(plot_path, dir)
    try mkdir(path) end

    # Limit to max disk and outside swarm
    dp = df[df[:InsideMaxDisk] & (df[:DistFromEdge].>=0),:]
    dp[:Spread] = dp[:ExpectedCascadeSize] - dp[:Detections]
    dp[:RelativeSpread] = dp[:Spread] ./ dp[:Detections]

    h = plot(dp, x=:DistFromEdge, y=:ExpectedCascadeSize,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Distance from boundary (body length)"),
        Guide.ylabel("Expected cascade size"),
        Guide.colorkey("State"))

    draw(PDF(joinpath(dir, "cascade_dist.pdf"), 6inch, 4inch), h)

    h = plot(dp, x=:DistFromEdge, y=:Spread,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Distance from boundary (body length)"),
        Guide.ylabel("Expected Spread"),
        Guide.colorkey("State"))

    draw(PDF(joinpath(dir, "spread_dist.pdf"), 6inch, 4inch), h)

    h = plot(dp, x=:DistFromEdge, y=:RelativeSpread,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Distance from boundary (body length)"),
        Guide.ylabel("Expected relative spread"),
        Guide.colorkey("State"))

    draw(PDF(joinpath(dir, "rel_spread_dist.pdf"), 6inch, 4inch), h)
end


function plot_spread_angle_dist(df::AbstractDataFrame; dir="spread")

    dir = joinpath(plot_path, dir)
    try mkdir(path) end

    # Limit to polarized, max disk and outside swarm
    idx = (df[:SwarmState].==:Polarized) & df[:InsideMaxDisk] & (df[:DistFromEdge].>=0)

    dp = df[idx,:]
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
        e = 0:5:100
        _, bin = hist(d[:ExpectedCascadeSize], e)
        DataFrame(ExpectedCascadeSize=(e[2:end]+e[1:end-1])/2, Density=bin ./ maximum(bin))
    end

    p = plot(dq, x=:BinRelativeDir, y=:ExpectedCascadeSize,
        color=:Density, Geom.rectbin,
        Coord.cartesian(xmin=-180, xmax=180),
        Scale.x_continuous(labels=x->@sprintf("%dº", x)),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Angle relative to group direction"),
        Guide.ylabel("Expected cascade size"),
        Guide.title("Polarized schools only"))

    draw(PDF(joinpath(dir, "spread_angle_dist.pdf"), 6inch, 4inch), p)
end

function plot_spread_rel_angle_dist(df::AbstractDataFrame; dir="spread")

    dir = joinpath(plot_path, dir)
    try mkdir(path) end

    # Limit to polarized, max disk and outside swarm
    idx = (df[:SwarmState].==:Polarized) & df[:InsideMaxDisk] & (df[:DistFromEdge].>=0)

    dp = df[idx,:]
    dp[:RelativeDir] = rad2deg(mod(dp[:RelativeDir] + 3pi, 2pi) - pi)
    dp[:BinRelativeDir] = convert(DataVector{Float64}, zeros(size(dp, 1)))
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])
    dp[dp[:Detections].==0,:Detections] = NA
    ey = -180:10:180
    for i in 1:size(dp, 1)
        for j in 2:length(ey)
            if dp[i,:RelativeDir] < ey[j]
                dp[i,:BinRelativeDir] = (ey[j-1] + ey[j]) / 2
                break
            end
        end
    end

    dq = by(dp, :BinRelativeDir) do d
        e = 0:0.05:1
        _, bin = hist((d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections], e)
        DataFrame(RelativeSpread=(e[2:end]+e[1:end-1])/2, Density=bin ./ maximum(bin))
    end

    p = plot(dq, x=:BinRelativeDir, y=:RelativeSpread,
        color=:Density, Geom.rectbin,
        Coord.cartesian(xmin=-180, xmax=180),
        Scale.x_continuous(labels=x->@sprintf("%dº", x)),
        Guide.xticks(ticks=collect(-180:30:180)),
        Guide.xlabel("Angle relative to group direction"),
        Guide.ylabel("Expected relative spread"),
        Guide.title("Polarized schools only"))

    draw(PDF(joinpath(dir, "rel_spread_angle_dist.pdf"), 6inch, 4inch), p)
end


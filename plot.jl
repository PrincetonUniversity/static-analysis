# Plotting related functions

using Gadfly
using Compose
using Colors

# export stdplot, plotalpha, netplot, plotscale, plotsize

function stdplot(df::AbstractDataFrame)
    # compute mean and SEM
    dfb = by(df[df[:SwarmState].==:Polarized,:], [:Kind, :BinFromBack]) do d
        m, s = mean(d[:Info]), sem(d[:Info])
        DataFrame(Dist=d[1,:BinDistFromBack], Mean=m, Min=m-s, Max=m+s)
    end

    # standardize units
    dfb = by(dfb, :Kind) do d
        m, s = mean(d[:Mean]), std(d[:Mean])
        DataFrame(Dist=d[:Dist], Mean=(d[:Mean]-m)./s, Min=(d[:Min]-m)./s, Max=(d[:Max]-m)./s)
    end

    # plot
    p1 = plot(unique(dfb), x=:Dist, y=:Mean, ymin=:Min, ymax=:Max, color=:Kind, Geom.point, Geom.errorbar,
        Guide.xlabel("Normalized distance from back to front (polarized groups only)"),
        Guide.ylabel("Standardized units"))

    # compute mean and SEM
    dfe = by(df[df[:SwarmState].==:Polarized,:], [:Kind, :BinFromEdge]) do d
        m, s = mean(d[:Info]), sem(d[:Info])
        d[1,:BinFromEdge] > 0 || return DataFrame()
        DataFrame(Dist=d[1,:BinDistFromEdge], Mean=m, Min=m-s, Max=m+s)
    end


    # standardize units
    dfe = by(dfe, :Kind) do d
        m, s = mean(d[:Mean]), std(d[:Mean])
        DataFrame(Dist=d[:Dist], Mean=(d[:Mean]-m)./s, Min=(d[:Min]-m)./s, Max=(d[:Max]-m)./s)
    end

    # plot
    p2 = plot(dfe, x=:Dist, y=:Mean, ymin=:Min, ymax=:Max, color=:Kind, Geom.point, Geom.errorbar,
        # Guide.xlabel("Normalized distance from edge to center"),
        Coord.cartesian(xmin=0, xmax=maximum(dfe[:Dist])+minimum(dfe[:Dist])),
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Standardized units"))

    p = vstack(p1, p2)

    # draw(PDF("info.pdf", 6inch, 8inch), p)

    return p
end

function bounds(p::Vector{State})
    xmin, xmax = extrema([v.pos.x for v in p])
    ymin, ymax = extrema([v.pos.y for v in p])
    (xmin - 1, xmax + 1, ymin - 1, ymax + 1)
end

"netplot plots a weighted directed network of oriented ellipses."
netplot(p::Vector{State}, w::Matrix{Float64}) = netplot(p, w, bounds(p))

function netplot(p::Vector{State}, w::Matrix{Float64}, bounds)
    xmin, xmax, ymin, ymax = bounds
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c1 = context()
    for i in eachindex(p), j in eachindex(p)
        w[i,j] > 0.01 || continue
        c1 = compose(c1, (context(), line([(p[i].pos.x, p[i].pos.y), (p[j].pos.x, p[j].pos.y)]), stroke(RGBA(i<j, 0, i>j, w[i,j]))))
    end
    c2 = context()
    for i in eachindex(p)
        ctx = context(rotation=Rotation(angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
        c2 = compose(c2, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2)))
        # c2 = compose(c2, myellipse(p[i].pos, Vec2(1, 0.125), 0.8, angle(p[i].vel)))
    end
    compose(context(units=box), c2, c1)
end

function myellipse(pos::Vec2, size::Vec2, offset::Real, θ::Real)
    N = 51
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    points = Vector{NTuple{2,Float64}}(N)
    for i in 1:N
        ϕ = 2(i-1)π / (N+1)
        x = R * [(cos(ϕ)-offset)*size.x/2, sin(ϕ)*size.y/2]
        points[i] = (pos.x + x[1], pos.y + x[2])
    end
    polygon(points)
end

colorplot(p::Vector{State}, by::Vector{Float64}) = colorplot(p, by, bounds(p))
function colorplot(p::Vector{State}, by::Vector{Float64}, bounds)
    xmin, xmax, ymin, ymax = bounds
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context()
    for i in eachindex(p)
        c = compose(c, (context(), myellipse(p[i].pos, Vec2(1, 0.125), 0.8, p[i].dir), fill(RGB(1, 1-by[i], 0))))
    end
    compose(context(units=box), c)
end

function imagedetect(p::Vector{State}, detect::Matrix{Int})
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

function plotdetect(p::Vector{State}, detect::Matrix{Int})
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

"netvid make a video of the dynamic network of interaction."
function netvid(p::Matrix{State}, mask::BitArray{2}, si::Array{Float64, 3})
    xmin, xmax, ymin, ymax = bounds(p[mask])
    r = (ymax - ymin) / (xmax - xmin)
    dir = "plots/vid/net"
    try rm(dir, recursive=true) end
    try mkdir(dir) end
    K = size(p, 2)
    for k in 1:K
        @printf("\r%3d%%", 100(k-1) / K)
        m = mask[:,k]
        n = netplot(p[m,k], si[m,m,k], (xmin, xmax, ymin, ymax))
        n = compose(context(), n, (context(), rectangle(), fill(colorant"white")))
        file = @sprintf("%08d.png", k)
        # draw(PNG(joinpath(dir, file), 8inch, r*8inch), n)
        draw(PNG(joinpath(dir, file), 1440px, 800px), n)
    end
    println("\r100%")
end

function colorvidperso(p::Matrix{State}, mask::BitArray{2}, pi::Array{Float64, 2})
    xmin, xmax, ymin, ymax = bounds(p[mask])
    r = (ymax - ymin) / (xmax - xmin)
    dir = "plots/vid/pi"
    try rm(dir, recursive=true) end
    try mkdir(dir) end
    K = size(p, 2)
    for k in 1:K
        @printf("\r%3d%%", 100(k-1) / K)
        m = mask[:,k]
        n = colorplot(p[m,k], pi[m,k], (xmin, xmax, ymin, ymax))
        n = compose(context(), n, (context(), rectangle(), fill(colorant"black")))
        file = @sprintf("%08d.png", k)
        # draw(PNG(joinpath(dir, file), 8inch, r*8inch), n)
        draw(PNG(joinpath(dir, file), 1440px, 800px), n)
    end
    println("\r100%")
end

function colorvidsocial(p::Matrix{State}, mask::BitArray{2}, si::Array{Float64, 3})
    xmin, xmax, ymin, ymax = bounds(p[mask])
    r = (ymax - ymin) / (xmax - xmin)
    dir = "plots/vid/si"
    try rm(dir, recursive=true) end
    try mkdir(dir) end
    K = size(p, 2)
    for k in 1:K
        @printf("\r%3d%%", 100(k-1) / K)
        m = mask[:,k]
        n = colorplot(p[m,k], 50000 .* cc(si[m,m,k]), (xmin, xmax, ymin, ymax))
        n = compose(context(), n, (context(), rectangle(), fill(colorant"black")))
        file = @sprintf("%08d.png", k)
        # draw(PNG(joinpath(dir, file), 8inch, r*8inch), n)
        draw(PNG(joinpath(dir, file), 1440px, 800px), n)
    end
    println("\r100%")
end


"plotalpha plots the precomputed α-shape of a group of particles."
function plotalpha(outer, inner, degen, solo, p::Vector{State})
    layers = Layer[]
    for u in outer
        v = [u; u[1]]
        x = [q.pos.x for q in p[v]]
        y = [q.pos.y for q in p[v]]
        append!(layers, layer(x=x, y=y, Geom.point, Geom.path, Theme(default_color=colorant"red")))
    end
    for u in inner
        v = [u; u[1]]
        x = [q.pos.x for q in p[v]]
        y = [q.pos.y for q in p[v]]
        append!(layers, layer(x=x, y=y, Geom.point, Geom.path, Theme(default_color=colorant"orange")))
    end
    for u in degen
        v = [u; u[1]]
        x = [q.pos.x for q in p[v]]
        y = [q.pos.y for q in p[v]]
        append!(layers, layer(x=x, y=y, Geom.point, Geom.path, Theme(default_color=colorant"gray")))
    end
    x = [q.pos.x for q in p[solo]]
    y = [q.pos.y for q in p[solo]]
    append!(layers, layer(x=x, y=y, Geom.point, Theme(default_color=colorant"pink")))
    x = [q.pos.x for q in p]
    y = [q.pos.y for q in p]
    append!(layers, layer(x=x, y=y, Geom.point))
    plot(layers...)
end

function plotmap(df::AbstractDataFrame)
    df2 = by(df[df[:SwarmState].==:Polarized,:], [:Kind, :BinFromBack, :BinFromEdge]) do d
        # if nrow(d) < 100 return DataFrame() end
        if d[1,:BinFromEdge] == 0 return DataFrame() end
        DataFrame(BinDistFromBack=d[1,:BinDistFromBack], BinDistFromEdge=d[1,:BinDistFromEdge], Info=median(d[:Info]), Count=nrow(d))
    end
    p = plot(df2[df2[:Kind].==:Personal,:], x=:BinDistFromBack, y=:BinDistFromEdge, color=:Info, Geom.rectbin,
        Guide.xlabel("Normalized distance from back to front"),
        Guide.ylabel("Distance from edge (body length)"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("External\nvisual field"))
    draw(PDF("heatmap_personal.pdf", 5.5inch, 4inch), p)
    p = plot(df2[df2[:Kind].==:Social,:], x=:BinDistFromBack, y=:BinDistFromEdge, color=:Info, Geom.rectbin,
        Guide.xlabel("Normalized distance from back to front"),
        Guide.ylabel("Distance from edge (body length)"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("Social influence"))
    draw(PDF("heatmap_social.pdf", 5.5inch, 4inch), p)

    df3 = by(df[df[:SwarmState].==:Polarized,:], [:Kind, :BinFromBack, :BinFromEdgeRel]) do d
        # if nrow(d) < 100 return DataFrame() end
        if d[1,:BinFromEdge] == 0 return DataFrame() end
        DataFrame(BinDistFromBack=d[1,:BinDistFromBack], BinDistFromEdge=d[1,:BinDistFromEdgeRel], Info=median(d[:Info]), Count=nrow(d))
    end
    p = plot(df3[df3[:Kind].==:Personal,:], x=:BinDistFromBack, y=:BinDistFromEdge, color=:Info, Geom.rectbin,
        Guide.xlabel("Normalized distance from back to front"),
        Guide.ylabel("Distance from edge (body length)"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("External\nvisual field"))
    draw(PDF("heatmap_personal_rel.pdf", 5.5inch, 4inch), p)
    p = plot(df3[df3[:Kind].==:Social,:], x=:BinDistFromBack, y=:BinDistFromEdge, color=:Info, Geom.rectbin,
        Guide.xlabel("Normalized distance from back to front"),
        Guide.ylabel("Distance from edge (body length)"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("Social influence"))
    draw(PDF("heatmap_social_rel.pdf", 5.5inch, 4inch), p)

    p = plot(df2[df2[:Kind].==:Personal,:], x=:BinDistFromBack, y=:BinDistFromEdge, color=:Count, Geom.rectbin,
        Guide.xlabel("Normalized distance from back to front"),
        Guide.ylabel("Distance from edge (body length)"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("Count"))
    draw(PDF("heatmap_personal_count.pdf", 5.5inch, 4inch), p)
    p = plot(df3[df3[:Kind].==:Personal,:], x=:BinDistFromBack, y=:BinDistFromEdge, color=:Count, Geom.rectbin,
        Guide.xlabel("Normalized distance from back to front"),
        Guide.ylabel("Distance from edge (body length)"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("Count"))
    draw(PDF("heatmap_social_count.pdf", 5.5inch, 4inch), p)
end

function plotscale(df::AbstractDataFrame)
    # df = readtable("data/scale/df.csv")
    # df[:Kind] = convert(DataArrays.DataArray{Symbol,1}, df[:Kind])
    # df[:SwarmState] = convert(DataArrays.DataArray{Symbol,1}, df[:SwarmState])
    ## df[:Scale] = convert(DataArrays.DataArray{Symbol,1}, map(symbol, df[:Scale]))

    # compute mean and SEM
    dfb = by(df[df[:SwarmState].==:Polarized,:], [:Kind, :BinFromBack, :Scale]) do d
        m, s = mean(d[:Info]), sem(d[:Info])
        DataFrame(Dist=d[1,:BinDistFromBack], Mean=m, Min=m-s, Max=m+s, Scale=d[1,:Scale])
    end

    # standardize units
    dfb = by(dfb, [:Kind, :Scale]) do d
        m, s = mean(d[:Mean]), std(d[:Mean])
        DataFrame(Dist=d[:Dist], Mean=(d[:Mean]-m)./s, Min=(d[:Min]-m)./s, Max=(d[:Max]-m)./s, Scale=d[:Scale])
    end


    # compute mean and SEM
    dfe = by(df, [:Kind, :BinFromEdge, :Scale]) do d
        m, s = mean(d[:Info]), sem(d[:Info])
        d[1,:BinFromEdge] > 0 || return DataFrame()
        DataFrame(Dist=d[1,:BinDistFromEdge], Mean=m, Min=m-s, Max=m+s, Scale=d[1,:Scale])
    end

    # standardize units
    dfe = by(dfe, [:Kind, :Scale]) do d
        m, s = mean(d[:Mean]), std(d[:Mean])
        DataFrame(Dist=d[:Dist], Mean=(d[:Mean]-m)./s, Min=(d[:Min]-m)./s, Max=(d[:Max]-m)./s, Scale=d[:Scale])
    end

    x0, x1 = extrema((linspace(0, 1, 12)[1:end-1] + linspace(0, 1, 12)[2:end])/2)
    dfe[:NormDist] = 0.0
    for scale in unique(dfe[:Scale])
        id = dfe[:Scale].==scale
        min, max = extrema(dfe[id,:Dist])
        dfe[id,:NormDist] = x0 + (dfe[id,:Dist] - min) * (x1 - x0) / (max - min)
    end

    p1 = plot(dfb[dfb[:Kind].==:Personal,:], x=:Dist, y=:Mean, ymin=:Min, ymax=:Max, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from back"),
        Guide.ylabel("Standardized units"),
        Guide.title("External visual field"),
        Scale.color_log2(minvalue=0.5, maxvalue=2.0))
    p2 = plot(dfe[dfe[:Kind].==:Personal,:], x=:Dist, y=:Mean, ymin=:Min, ymax=:Max, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Standardized units"),
        Guide.title("External visual field"),
        Scale.color_log2(minvalue=0.5, maxvalue=2.0))
    p3 = plot(dfe[dfe[:Kind].==:Personal,:], x=:NormDist, y=:Mean, ymin=:Min, ymax=:Max, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from edge"),
        Guide.ylabel("Standardized units"),
        Guide.title("External visual field"),
        Scale.color_log2(minvalue=0.5, maxvalue=2.0))
    p = vstack(p1, p2, p3)

    draw(PDF("scale_personal.pdf", 6inch, 12inch), p)

    p1 = plot(dfb[dfb[:Kind].==:Social,:], x=:Dist, y=:Mean, ymin=:Min, ymax=:Max, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from back"),
        Guide.ylabel("Standardized units"),
        Guide.title("Social influence"),
        Scale.color_log2(minvalue=0.5, maxvalue=2.0))
    p2 = plot(dfe[dfe[:Kind].==:Social,:], x=:Dist, y=:Mean, ymin=:Min, ymax=:Max, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Standardized units"),
        Guide.title("Social influence"),
        Scale.color_log2(minvalue=0.5, maxvalue=2.0))
    p3 = plot(dfe[dfe[:Kind].==:Social,:], x=:NormDist, y=:Mean, ymin=:Min, ymax=:Max, color=:Scale, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from edge"),
        Guide.ylabel("Standardized units"),
        Guide.title("Social influence"),
        Scale.color_log2(minvalue=0.5, maxvalue=2.0))
    p = vstack(p1, p2, p3)

    draw(PDF("scale_social.pdf", 6inch, 12inch), p)
end

function plotsize(mode)
    dfb = readtable("size/dfb_$mode.csv")
    dfb[:Kind] = convert(DataArrays.DataArray{Symbol,1}, dfb[:Kind])
    dfe = readtable("size/dfe_$mode.csv")
    dfe[:Kind] = convert(DataArrays.DataArray{Symbol,1}, dfe[:Kind])

    x0, x1 = extrema((linspace(0, 1, 12)[1:end-1] + linspace(0, 1, 12)[2:end])/2)
    dfe[:NormDistFromEdge] = 0.0
    for size in unique(dfe[:Size])
        id = dfe[:Size].==size
        min, max = extrema(dfe[id,:DistFromEdge])
        dfe[id,:NormDistFromEdge] = x0 + (dfe[id,:DistFromEdge] - min) * (x1 - x0) / (max - min)
    end

    p1 = plot(dfb[dfb[:Kind].==:Personal,:], x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from back"),
        Guide.ylabel("Standardized units"),
        Guide.title("External visual field ($mode)"),
        Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
    p2 = plot(dfe[dfe[:Kind].==:Personal,:], x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Standardized units"),
        Guide.title("External visual field ($mode)"),
        Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
    p3 = plot(dfe[dfe[:Kind].==:Personal,:], x=:NormDistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from edge"),
        Guide.ylabel("Standardized units"),
        Guide.title("External visual field ($mode)"),
        Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
    p = vstack(p1, p2, p3)

    draw(PDF("size_personal_$(mode).pdf", 6inch, 12inch), p)

    p1 = plot(dfb[dfb[:Kind].==:Social,:], x=:DistFromBack, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from back"),
        Guide.ylabel("Standardized units"),
        Guide.title("Social influence ($mode)"),
        Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
    p2 = plot(dfe[dfe[:Kind].==:Social,:], x=:DistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Distance from edge (body length)"),
        Guide.ylabel("Standardized units"),
        Guide.title("Social influence ($mode)"),
        Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
    p3 = plot(dfe[dfe[:Kind].==:Social,:], x=:NormDistFromEdge, y=:MeanInfo, ymin=:LowInfoSE, ymax=:HighInfoSE, color=:Size, Geom.point, Geom.line, Geom.errorbar,
        Guide.xlabel("Normalized distance from edge"),
        Guide.ylabel("Standardized units"),
        Guide.title("Social influence ($mode)"),
        Scale.color_continuous(minvalue=5.0, maxvalue=25.0))
    p = vstack(p1, p2, p3)

    draw(PDF("size_social_$(mode).pdf", 6inch, 12inch), p)
end

function plot_detect_radius(df::AbstractDataFrame)
    # Limit to max disk and outside swarm
    dp = df[df[:InsideMaxDisk] & (df[:DistFromEdge].>=0),:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])

    return plot(dp, x=:DistFromEdge, y=:Detections,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Distance from boundary (body length)"),
        Guide.ylabel("Number of overlapping fields of view"),
        Guide.colorkey("State"))
end

function plot_detect_polar(df::AbstractDataFrame)
    # Limit to polarized, max disk and outside swarm
    dp = df[(df[:SwarmState].==:Polarized) & df[:InsideMaxDisk] & (df[:DistFromEdge].>=0),:]
    dp[:Detections] = convert(DataVector{Float64}, dp[:Detections])
    dp[:RelativeDir] = mod(dp[:RelativeDir] + 3pi, 2pi) - pi

    return plot(dp, x=:RelativeDir, y=:Detections,
        color=:SwarmState, Geom.line, Stat.binmean(n=50),
        Guide.xlabel("Angle relative to group direction"),
        Guide.ylabel("Number of overlapping fields of view"),
        Guide.title("Polarized groups only"),
        Guide.colorkey("State"))
end

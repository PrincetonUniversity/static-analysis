# Compute personal and social information from ellipswarm/static.

@step [get_net_file] function get_personal_info(p::Project)
    h5read(p.step[get_net_file], "personal")
end

@step [get_net_file] function get_social_info(p::Project)
    h5read(p.step[get_net_file], "social")
end

@step [get_pos, get_personal_info, get_social_info, run_dist_from_edge] function run_info(p::Project)
    pos, mask = p.step[get_pos]
    α = p.conf[:α]
    pi = p.step[get_personal_info]
    si = p.step[get_social_info]

    N = size(pos, 1) # max swarm size
    K = size(pos, 2) # replicates

    # distance from edge
    dem, msk = p.step[run_dist_from_edge]
    de_max = maximum(dem)

    # distances, order parameters, dataframe
    cols = [
        (:DistFromBack,       Float64),
        (:DistFromEdge,       Float64),
        (:DistFromEdgeRel,    Float64),
        (:SwarmID,            Int),
        (:Polarization,       Float64),
        (:Rotation,           Float64),
        (:State,              ASCIIString),
        (:Info,               Float64),
        (:Kind,               ASCIIString),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k] & msk[:,k]

        # distances
        db = dist_from_back(pos[nz,k])
        de = dem[nz,k]
        de_rel = de ./ maximum(de)

        # order parameters and state
        op, or = order_parameters(pos[nz,k])
        if op > 0.65 && or < 0.35
            state = "Polarized"
        elseif op < 0.35 && or > 0.65
            state = "Milling"
        elseif op < 0.35 && or < 0.35
            state = "Swarming"
        else
            state = "Transitioning"
        end

        # social influence
        msi = cc(si[nz,nz,k])

        # append to dataframe
        row = DataFrame(
            DistFromBack       = db,
            DistFromEdge       = de,
            DistFromEdgeRel    = de_rel,
            SwarmID            = k,
            Polarization       = op,
            Rotation           = or,
            State              = state,
            Info               = pi[nz,k],
            Kind               = "Personal",
        )
        append!(df, row)
        row[:Info] = msi
        row[:Kind] = "Social"
        append!(df, row)
    end
    pool!(df, [:State, :Kind])
    println("\r100%")

    df
end


# Plots
# -----

@step [run_info] function plot_evf_si_std(p::Project)
    df = copy(p.step[run_info])
    theme = copy(p.conf[:theme])
    theme[:major_label_font_size] = 10pt * 4/3
    theme[:minor_label_font_size] = 9pt * 4/3
    theme[:key_title_font_size] = 0pt
    theme[:key_label_font_size] = 9pt * 4/3

    # polarized groups only as in Rosenthal et al. (2015) PNAS
    df = df[df[:State].=="Polarized",:]

    function plot_std(col, upper, label)

        # bin distances
        nbins = 11
        bin_edge = linspace(0, upper, nbins + 1)
        for row in eachrow(df)
            for j in 2:length(bin_edge)
                if row[col] < bin_edge[j]
                    row[col] = (bin_edge[j-1] + bin_edge[j]) / 2
                    break
                end
            end
        end
        sub = df[col] .< upper

        # compute mean and SEM
        dfe = by(df[sub,:], [:Kind, col]) do d
            m, s = mean(d[:Info]), sem(d[:Info])
            DataFrame(Mean=m, Min=m-s, Max=m+s)
        end

        # standardize units
        dfe = by(dfe, :Kind) do d
            m, s = mean(d[:Mean]), std(d[:Mean]) * 2
            DataFrame(Dist=d[col], Mean=(d[:Mean]-m)./s, Min=(d[:Min]-m)./s, Max=(d[:Max]-m)./s)
        end

        # plot
        color_social = colorant"#f12"
        color_personal = colorant"#08c"
        p = plot(
            layer(dfe[dfe[:Kind].=="Social",:],
                x=:Dist, ymin=:Min, ymax=:Max,
                Geom.errorbar, order=2,
                Theme(default_color=color_social)),
            layer(dfe[dfe[:Kind].=="Social",:],
                x=:Dist, y=:Mean,
                Geom.point, order=3,
                Theme(default_color=color_social)),
            layer(dfe[dfe[:Kind].=="Personal",:],
                x=:Dist, ymin=:Min, ymax=:Max,
                Geom.errorbar, order=0,
                Theme(default_color=color_personal)),
            layer(dfe[dfe[:Kind].=="Personal",:],
                x=:Dist, y=:Mean,
                Geom.point, order=1,
                Theme(default_color=color_personal)),
            Coord.cartesian(xmin=0, xmax=upper, ymin=-0.7, ymax=1.25),
            Guide.xlabel(label),
            Guide.ylabel("Standardized units"),
            Guide.manual_color_key("",
                ["social influence", "external visual field"],
                [color_social, color_personal]),
            Theme(key_position=:top; theme...))

        name = col == :DistFromBack ? "evf_si_std_back.pdf" : "evf_si_std_edge.pdf"
        draw(PDF(joinpath(plot_path, name), 4inch, 3inch), p)
    end

    print("  0%")
    max_dist_from_edge = 5
    plot_std(:DistFromEdge, max_dist_from_edge,
        "Distance from edge (body length)")
    print("\r 50%")
    max_dist_from_back = 1
    plot_std(:DistFromBack, max_dist_from_back,
        "Normalized distance from back to front")
    println("\r100%")
end


@step [get_pos, get_personal_info, get_social_info] function plot_info_illustration(p::Project)
    pos, mask = p.step[get_pos]
    pi = p.step[get_personal_info]
    si = p.step[get_social_info]
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    N = size(pos, 1) # max swarm size
    K = size(pos, 2) # replicates

    k = 1 # selected school ID

    m = mask[:,k]

    xmin, xmax, ymin, ymax = bounds(pos[m,k])
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)

    function plot_color(info, label, tag, scale, cmin, cmax)
        c = context()
        for i in eachindex(pos[m,k])
            p = pos[m,k]
            ctx = context(rotation=Rotation(-angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
            c = compose(c, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2), fill(cmap(info[i]))))
        end
        c = compose(context(units=box), c)
        p = plot(x=[100,101], y=[100,101], color=[0.,1], Geom.point,
            Guide.annotation(c),
            Guide.xlabel("x (body length)"),
            Guide.ylabel("y (body length)"),
            Guide.colorkey(label),
            Theme(; base_theme...),
            scale(minvalue=cmin, maxvalue=cmax, colormap=cmap),
            Coord.cartesian(raster=true, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fixed=true))
        path = joinpath(plot_path, string("info_illustration_", tag, ".pdf"))
        draw(PDF(path, 5inch, 4inch), p)
    end

    plot_color(pi[m,k], "external\nvisual field", "evf",
        Scale.color_continuous, 0, 1)
    rescale(x) = (log10(x) + 6) / 4
    plot_color(rescale(cc(si[m,m,k])), "social\ninfluence", "si",
        Scale.color_log10, 1e-6, 1e-2)
end


@step [run_info] function plot_info_dist(p::Project)
    df = copy(p.step[run_info])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to outside school boundaries
    df = df[df[:DistFromEdge] .>= 0,:]

    # Rescale personal info from [0,1] to [0,360]
    df[df[:Kind] .== "Personal", :Info] *= 360

    function bin_dist(dist, edges)
        for row in eachrow(df)
            for j in 2:length(edges)
                if row[dist] < edges[j]
                    row[dist] = (edges[j-1] + edges[j]) / 2
                    break
                end
            end
        end
        df[df[dist] .< last(edges),:]
    end

    function plot_dist(df, dist, info, tag, xlabel, xmax, ylabel, ymin, ymax, edges, ticks, scale)
        dfa = by(df[df[:Kind] .== info, :], dist) do d
            _, bin = hist(d[:Info], edges)
            DataFrame(Info=(edges[2:end]+edges[1:end-1])/2, Density=bin ./ maximum(bin))
        end

        dfm = by(df[df[:Kind] .== info, :], dist) do d
            DataFrame(Info=median(d[:Info]))
        end

        p = plot(
            layer(dfa, x=dist, y=:Info, color=:Density, Geom.rectbin),
            layer(dfm, x=dist, y=:Info, Geom.line, order=1,
                Theme(default_color=colorant"black")),
            Coord.cartesian(xmin=0, xmax=xmax, ymin=ymin, ymax=ymax),
            scale,
            Scale.color_continuous(colormap=cmap),
            ticks,
            Guide.xlabel(xlabel),
            Guide.ylabel(ylabel),
            Theme(; base_theme...))

        name = string("info_dist_", tag, ".pdf")
        draw(PDF(joinpath(plot_path, name), 6inch, 4inch), p)
    end

    print("  0%")

    dfe = bin_dist(:DistFromEdge, 0:0.5:5)
    print("\r 40%")

    plot_dist(dfe, :DistFromEdge, "Personal", "evf_edge",
        "Distance from edge (body length)", 5,
        "External visual field", 0, 360, 0:5:360,
        Guide.yticks(ticks=collect(0:45:360)),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)))
    print("\r 45%")

    plot_dist(dfe, :DistFromEdge, "Social", "si_edge",
        "Distance from edge (body length)", 5,
        "Social influence", -6, -2, logspace(-6, -2, 64),
        Guide.yticks,
        Scale.y_log10)
    print("\r 50%")

    dfb = bin_dist(:DistFromBack, 0:0.025:1)
    print("\r 90%")

    plot_dist(dfb, :DistFromBack, "Personal", "evf_back",
        "Normalized distance from back to front", 1,
        "External visual field", 0, 360, 0:5:360,
        Guide.yticks(ticks=collect(0:45:360)),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)))
    print("\r 95%")

    plot_dist(dfb, :DistFromBack, "Social", "si_back",
        "Normalized distance from back to front", 1,
        "Social influence", -6, -2, logspace(-6, -2, 64),
        Guide.yticks,
        Scale.y_log10)
    println("\r100%")
end


@step [plot_evf_si_std, plot_info_illustration, plot_info_dist] function plot_chap2(p::Project)
    println("All done")
end


@step [get_pos] function plot_example_school(p::Project)
    pos, mask = p.step[get_pos]
    k = 1
    q = pos[mask[:,k],k];
    draw(PDF(joinpath(plot_path, "example.pdf"), 6inch, 6inch), simpleplot(q))
end


function plot_weights(p::Project...)
    cols = [
        (:Weight,       Float64),
        (:Project,      ASCIIString),
        (:SwarmID,      Int),
        (:Polarization, Float64),
        (:Rotation,     Float64),
        (:State,        ASCIIString),
    ]
    df = DataFrame([x[2]::DataType for x in cols], [x[1]::Symbol for x in cols], 0)
    pool!(df, [:Project, :State])
    for q in p
        println(q.uid, ":")
        pos, mask = get_pos(q)
        si = get_social_info(q)

        K = size(pos, 2) # replicates
        for k in 1:K
            @printf("\r%3d%%", 100(k-1) / K)
            nz = mask[:,k]
            op, or = order_parameters(pos[nz,k])
            if op > 0.65 && or < 0.35
                state = "Polarized"
            elseif op < 0.35 && or > 0.65
                state = "Milling"
            elseif op < 0.35 && or < 0.35
                state = "Swarming"
            else
                state = "Transitioning"
            end
            net = si[nz,nz,k]
            row = DataFrame(
                Weight          = net[net.>0],
                Project         = q.uid,
                SwarmID         = k,
                Polarization    = op,
                Rotation        = or,
                State           = state,
            )
            append!(df, row)
        end
        println("\r100%")
    end

    h = plot(df, x=:Weight, ygroup=:State, color=:Project,
        Geom.subplot_grid(Geom.histogram),
        Scale.x_log10,
        Coordinate(xmax=0),
        Theme(; p[1].conf[:theme]...))

    draw(PDF("weights.pdf", 6inch, 8inch), h)

    df
end





# Helpers
# -------

"cc computes the local weighted directed clustering coefficient."
function cc(w::Matrix{Float64})
    # assuming w is square weight matrix
    n = size(w, 1)
    a = w .!= 0 # adjacency matrix
    CC = zeros(n)
    for i=1:n
        dg, db, C = 0.0, 0.0, 0.0
        for j=1:n
            for k=1:n
                C += (w[i,j] + w[j,i]) * (w[i,j] + w[k,i]) * (w[j,k] + w[k,j])
            end
            j != i || continue
            dg += a[i,j] + a[j,i]
            db += a[i,j] * a[j,i]
        end
        CC[i] = C == 0 ? 0 : C / 2(dg * (dg - 1) - 2db)
    end
    CC
end

"dist_from_back computes the distance of each particle to the back of the swarm."
function dist_from_back(p::Vector{State})
    m = mean(p)
    db = similar(p, Float64)
    for i in eachindex(p)
        @inbounds db[i] = dot(p[i].pos - m.pos, p[i].vel)
    end
    min, max = extrema(db)
    (db - min) ./ (max - min)
end


"order_parameters computes the polarization and rotation of a swarm."
function order_parameters(p::Vector{State})
    m = mean(p)
    up, ur = Vec2(0, 0), 0.0
    @inbounds for i in eachindex(p)
        u = unit(p[i].vel)
        up += u
        ur += cross(u, unit(p[i].pos - m.pos))
    end
    (norm(up) / length(p), abs(ur) / length(p))
end


function run_order_parameters(file::AbstractString)
    p, mask = h5read_particles(file)
    K = size(p, 2) # replicates
    op = zeros(K)
    or = zeros(K)
    for k=1:K
        @printf("\r%3d%%", 100(k-1) / K)
        nz = mask[:,k]
        op[k], or[k] = order_parameters(p[nz,k])
    end
    println("\r100%")
    op, or
end


function gen_dist_alpha(file::AbstractString)
    for r in [1, 2, 5, 10, 20, Inf]
        path = joinpath(data_path, string("dist_from_edge_", r, ".h5"))
        println(path, ":")
        de, mask = dist_from_edge(file, α=-1/r)
        try rm(path) end
        h5write(path, "DistFromEdge", de)
        h5write(path, "Mask", convert(Matrix{Int8}, mask))
    end
end


function plot_max_si(file::AbstractString; threshold=1e-3)
    p, mask = h5read_particles(file)
    si = h5read(file, "social")

    N = size(p, 1) # max swarm size
    K = size(p, 2) # replicates

    dem = h5read(joinpath(data_path, "dist_from_edge_5.0.h5"), "DistFromEdge")
    msk = h5read(joinpath(data_path, "dist_from_edge_5.0.h5"), "Mask") .!= 0

    for k=1:K
        nz = mask[:,k] & msk[:,k]
        msi = cc(si[nz,nz,k])
        de = dem[nz,k]
        # h = colorplot(p[nz,k], 0.25(msi .> threshold) + 0.75((msi .> threshold) & (de .> 0.03125/2) & (de .< 0.5)), (-15,15,-15,15))
        h = colorplot(p[nz,k], 0.25((msi .> 2e-4) & (msi .< 1e-3)) + 0.75((msi .> 2e-4) & (msi .< 1e-3) & (de .> 0.03125/2) & (de .<= 1.5*0.03125)), (-15,15,-15,15))
        draw(PNG(@sprintf("max_si/%06d.png", k), 900px, 900px), h)
    end
end


# Other plots
# -----------


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
        ctx = context(rotation=Rotation(-angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
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

colorplot{T<:Real}(p::Vector{State}, by::AbstractVector{T}) = colorplot(p, by, bounds(p))
function colorplot{T<:Real}(p::Vector{State}, by::AbstractVector{T}, bounds)
    xmin, xmax, ymin, ymax = bounds
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context()
    for i in eachindex(p)
        ctx = context(rotation=Rotation(-angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
        c = compose(c, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2), fill(RGB(by[i], by[i], 0.5-by[i]/2))))
        # c = compose(c, (context(), myellipse(p[i].pos, Vec2(1, 0.125), 0.8, p[i].dir), fill(RGB(1, 1-by[i], 0))))
    end
    compose(context(units=box), c)
end


simpleplot(p::Vector{State}) = simpleplot(p, bounds(p))
function simpleplot(p::Vector{State}, bounds)
    xmin, xmax, ymin, ymax = bounds
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context()
    for i in eachindex(p)
        ctx = context(rotation=Rotation(-angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
        c = compose(c, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2)))
        # c = compose(c, (context(), myellipse(p[i].pos, Vec2(1, 0.125), 0.8, p[i].dir)))
    end
    compose(context(units=box), c)
end


function distance_from_edge_plot(p::Vector{State}, de::Matrix{Float64}, grid::Matrix{Vec2})
    plot(x=[z.x for z in grid[:,1]], y=[z.y for z in grid[1,:]], z=de, Geom.contour(levels=-5:0.25:5),
        Guide.annotation(simpleplot(p, (-15, 15, -15, 15))),
        Guide.xlabel("x (body length)"),
        Guide.ylabel("y (body length)"),
        Coord.cartesian(xmin=-15, xmax=15, ymin=-15, ymax=15, fixed=true)
    )
end

function distance_from_edge_plot2(p::Vector{State};
                  minvalue=-5, maxvalue=5, title="Distance from edge (body length)",
                  colormap=x -> (x = minvalue + x * (maxvalue - minvalue); RGB(x < 0 ? 0 : 1, mod(x, 0.5), x < 0 ? 0 : mod(x, 0.5))))
    xmin, xmax = -15, 15
    ymin, ymax = -15, 15
    grid = [Vec2(x,y) for x in xmin:0.1:xmax, y in ymin:0.1:ymax]
    nx, ny = size(grid)
    x = [z.x for z in grid[:]]
    y = [z.y for z in grid[:]]
    de = dist_from_edge(0, p, grid)
    box = UnitBox(xmin, ymax, xmax - xmin, ymin - ymax)
    c = context()
    for i in eachindex(p)
        ctx = context(rotation=Rotation(-angle(p[i].vel), (p[i].pos.x, p[i].pos.y)))
        c = compose(c, (ctx, ellipse(p[i].pos.x - 0.8/2, p[i].pos.y, 1/2, 0.125/2)))
    end
    plot(x=x, y=y, color=de, Geom.rectbin,
        Guide.annotation(compose(context(), c, fill(colorant"gold"))),
        Guide.xlabel("x (body length)"),
        Guide.ylabel("y (body length)"),
        Guide.title(title),
        Scale.color_continuous(minvalue=minvalue, maxvalue=maxvalue, colormap=colormap),
        Coord.cartesian(raster=true, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fixed=true))
end


"netvid make a video of the dynamic network of interaction."
function netvid(p::Matrix{State}, mask::BitArray{2}, si::Array{Float64, 3})
    xmin, xmax, ymin, ymax = bounds(p[mask])
    r = (ymax - ymin) / (xmax - xmin)
    dir = joinpath(plot_path, "vidnet")
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
    dir = joinpath(plot_path, "vidpi")
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
    dir = joinpath(plot_path, "vidsi")
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


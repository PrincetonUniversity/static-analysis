# Analyse effect of changing interindividual distances.


@step function run_rescaling_info(p::Project)
    rps = rescaled_projects()
    for rp in rps
        println(format("\n# Scale {:.2f}", rp.conf[:scale]))
        run_info(rp)
        rp.step[run_info][:Scale] = rp.conf[:scale]
    end
    vcat(DataFrame[rp.step[run_info] for rp in rps])
end


@step function run_rescaling_structure(p::Project)
    rps = rescaled_projects()
    for rp in rps
        println(format("\n# Scale {:.2f}", rp.conf[:scale]))
        run_structure(rp)
        rp.step[run_structure][:Scale] = rp.conf[:scale]
    end
    vcat(DataFrame[rp.step[run_structure] for rp in rps])
end


@step function run_rescaling_spread(p::Project)
    rps = rescaled_sub_projects()
    for rp in rps
        println(format("\n# Scale {:.2f} - Sensitivity {:.1f}", rp.conf[:scale], rp.conf[:sensitivity]))
        run_spread(rp)
        rp.step[run_spread][:Scale] = rp.conf[:scale]
        rp.step[run_spread][:Sensitivity] = rp.conf[:sensitivity]
    end
    vcat(DataFrame[rp.step[run_spread] for rp in rps])
end

@step function run_rescaling_spread_summary(p::Project)
    rps = rescaled_sub_projects()
    for rp in rps
        println(format("\n# Scale {:.2f} - Sensitivity {:.1f}", rp.conf[:scale], rp.conf[:sensitivity]))
        run_spread_summary(rp)
    end
    vcat(DataFrame[rp.step[run_spread_summary] for rp in rps])
end

function load_run_spread(p::Project)
    dir = joinpath(data_path, "Projects", p.uid)
    file = joinpath(dir, "run_spread.jld")
    JLD.load(file, "data")
end

# TODO: save only summary (not 500MB!)
@step function run_spread_summary(p::Project)
    print("  0%")
    df = load_run_spread(p)
    print("\r 10%")

    # # Convert relative dir to degrees
    # df[:RelativeDir] = rad2deg(mod(df[:RelativeDir] + 3pi, 2pi) - pi)
    #
    # print("\r 20%")
    #
    # df[:Detections] = convert(DataVector{Float64}, df[:Detections])
    #
    # print("\r 30%")
    #
    # ey = -180:15:180
    # for i in 1:size(df, 1)
    #     for j in 2:length(ey)
    #         if df[i,:RelativeDir] < ey[j]
    #             df[i,:RelativeDir] = (ey[j-1] + ey[j]) / 2
    #             break
    #         end
    #     end
    # end

    print("\r 50%")

    # dfm = by(df, :RelativeDir) do d
    #     rho = d[:ExpectedCascadeSize]
        rho = df[:ExpectedCascadeSize]
        m1, s1 = mean(rho[!isnan(rho)]), sem(rho[!isnan(rho)])
        rho -= (df[:Detections])
        m2, s2 = mean(rho[!isnan(rho)]), sem(rho[!isnan(rho)])
        rho ./= (df[:Detections])
        m3, s3 = mean(rho[!isnan(rho)]), sem(rho[!isnan(rho)])
        dfm = DataFrame(
            # State=d[:State],
            # MaxRadius=d[:MaxRadius],
            # DistFromEdge=d[:DistFromEdge],
            Scale=p.conf[:scale],
            Sensitivity=p.conf[:sensitivity],
            ExpectedCascadeSize=m1,
            SEMCascadeSize=s1,
            ExpectedSpread=m2,
            SEMSpread=s2,
            ExpectedRelativeSpread=m3,
            SEMRelativeSpread=s3,
        )
    # end
    println("\r100%")
    dfm
end

@step [run_rescaling_info] function plot_rescaling_info(p::Project)
    print("  0%")

    df = copy(p.step[run_rescaling_info])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized swarms
    df = df[df[:State].=="Polarized",:]

    # Rescale personal info from [0,1] to [0,360]
    df[df[:Kind] .== "Personal", :Info] *= 360

    # Normalized distance from edge
    df[:NormDistFromEdge] = df[:DistFromEdge] ./ df[:Scale]

    # Round scaling factors
    for scale in [1/2,1/√2,1,√2,2]
        df[df[:Scale].==scale,:Scale] = round(scale, 2)
    end

    function bin_dist(dist, edges)
        @inbounds for row in eachrow(df)
            for j in 2:length(edges)
                if row[dist] < edges[j]
                    row[dist] = (edges[j-1] + edges[j]) / 2
                    break
                end
            end
        end
        df[df[dist] .< last(edges),:]
    end

    function plot_dist(df, dist, info, tag, xlabel, xmax, ylabel, ymin, ymax, ticks, scale)
        # compute mean and SEM
        # sub = (df[:Kind] .== info) & (df[:Info] .!= Inf)
        sub = Bool[(row[:Kind] == info) & (row[:Info] != Inf) for row in eachrow(df)]
        dfm = by(df[sub,:], [:Scale, dist]) do d
            if info == "Social"
                m, s = median(d[:Info]), sem(log10(d[:Info]))
                DataFrame(Median=m, Min=m*10^-s, Max=m*10^s)
            else
                m, s = median(d[:Info]), sem(d[:Info])
                DataFrame(Median=m, Min=m-s, Max=m+s)
            end
        end

        h = plot(dfm, x=dist, y=:Median, ymin=:Min, ymax=:Max, color=:Scale,
            Geom.point, Geom.line, Geom.errorbar,
            # Scale.color_log2(colormap=cmap),
            Scale.color_discrete_manual(map(cmap, linspace(0, 0.9, 5))...),
            Coord.cartesian(xmin=0, xmax=xmax, ymin=ymin, ymax=ymax),
            scale,
            ticks,
            Guide.xlabel(xlabel),
            Guide.ylabel(ylabel),
            Theme(; base_theme...))

        name = string("rescaling_", tag, ".pdf")
        draw(PDF(joinpath(plot_path, name), 6inch, 4inch), h)
    end

    print("\r 10%")

    dfe = bin_dist(:DistFromEdge, linspace(0, 5, 11))
    print("\r 20%")

    plot_dist(dfe, :DistFromEdge, "Personal", "evf_edge",
        "Distance from edge (body length)", 5,
        "External visual field", 0, 360,
        Guide.yticks(ticks=collect(0:45:360)),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)))
    print("\r 30%")

    plot_dist(dfe, :DistFromEdge, "Social", "si_edge",
        "Distance from edge (body length)", 5,
        "Social influence", -7, -2,
        Guide.yticks,
        Scale.y_log10)
    print("\r 40%")

    dfb = bin_dist(:DistFromBack, linspace(0, 1, 11))
    print("\r 50%")

    plot_dist(dfb, :DistFromBack, "Personal", "evf_back",
        "Normalized distance from back to front", 1,
        "External visual field", 0, 360,
        Guide.yticks(ticks=collect(0:45:360)),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)))
    print("\r 60%")

    plot_dist(dfb, :DistFromBack, "Social", "si_back",
        "Normalized distance from back to front", 1,
        "Social influence", -7, -2,
        Guide.yticks,
        Scale.y_log10)
    print("\r 70%")

    dfn = bin_dist(:NormDistFromEdge, linspace(0, 5, 11))
    print("\r 80%")

    plot_dist(dfn, :NormDistFromEdge, "Personal", "evf_edge_norm",
        "Normalized distance from edge", 5,
        "External visual field", 0, 360,
        Guide.yticks(ticks=collect(0:45:360)),
        Scale.y_continuous(labels=x->@sprintf("%d°", x)))
    print("\r 90%")

    plot_dist(dfn, :NormDistFromEdge, "Social", "si_edge_norm",
        "Normalized distance from edge", 5,
        "Social influence", -7, -2,
        Guide.yticks,
        Scale.y_log10)
    println("\r100%")
end


@step [run_rescaling_structure] function plot_rescaling_structure(p::Project)
    print("  0%")

    df = copy(p.step[run_rescaling_structure])
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Filter neighbor angles
    df = df[isnan(df[:NeighborAngle]),:]

    # Convert angles from rad to deg
    df[:ComponentAngle] *= 180/pi

    # Round scaling factors
    for scale in [1/2,1/√2,1,√2,2]
        df[df[:Scale].==scale,:Scale] = round(scale, 2)
    end

    # Filter components with a well-defined orientation
    wdo = df[:ComponentAspectRatio] .< 1/√2

    function plot_col(df, col, edges, tag, ticks, scale, ylabel)
        dfm = by(df, [:WeakestLink, :Scale]) do d
            if col == :ComponentAngle
                pos, neg = d[d[col].>0,col], d[d[col].<0,col]
                mpos, spos = median(pos), sem(pos)
                mneg, sneg = median(neg), sem(neg)
                DataFrame(
                    MedianPos=mpos, MinPos=mpos-spos, MaxPos=mpos+spos,
                    MedianNeg=mneg, MinNeg=mneg-sneg, MaxNeg=mneg+sneg)
            else
                m, s = median(d[col]), sem(d[col])
                DataFrame(Median=median(d[col]), Min=m-s, Max=m+s)
            end
        end
        layers = Layer[]
        if col == :ComponentSize
            append!(layers, [
                layer(dfm, x=:WeakestLink, ymin=:Min, ymax=:Max, color=:Scale,
                    Geom.errorbar, order=0);
                layer(dfm, x=:WeakestLink, y=:Median, color=:Scale,
                    Geom.line, order=1);
                layer(dfm, x=:WeakestLink, y=:Median, color=:Scale,
                    Geom.point, order=2)
            ])
        elseif col == :ComponentAngle
            append!(layers, [
            layer(dfm, x=:WeakestLink, ymin=:MinPos, ymax=:MaxPos, color=:Scale,
                Geom.errorbar, order=0);
            layer(dfm, x=:WeakestLink, y=:MedianPos, color=:Scale,
                Geom.line, order=1);
            layer(dfm, x=:WeakestLink, y=:MedianPos, color=:Scale,
                Geom.point, order=2);
            layer(dfm, x=:WeakestLink, ymin=:MinNeg, ymax=:MaxNeg, color=:Scale,
                Geom.errorbar, order=0);
            layer(dfm, x=:WeakestLink, y=:MedianNeg, color=:Scale,
                Geom.line, order=1);
            layer(dfm, x=:WeakestLink, y=:MedianNeg, color=:Scale,
                Geom.point, order=2)
            ])
        end
        h = plot(
            layers...,
            Coord.cartesian(xmin=-2.5, xmax=-0.5, ymin=edges[1], ymax=edges[end]),
            scale,
            Scale.x_log10,
            Scale.color_discrete_manual(map(cmap, linspace(0, 0.9, 5))...),
            ticks,
            Guide.xlabel("Weakest link threshold"),
            Guide.ylabel(ylabel),
            Theme(; base_theme...))

        name = string("rescaling_", tag, ".pdf")
        draw(PDF(joinpath(plot_path, name), 6inch, 4inch), h)
    end

    print("\r 40%")

    plot_col(df, :ComponentSize, 0:5:160, "component_size",
        Scale.y_continuous,
        Guide.yticks,
        "Component size")
    print("\r 70%")

    plot_col(df[wdo,:], :ComponentAngle, -90:5:90, "component_angle",
        Scale.y_continuous(labels=x->@sprintf("%d°", x)),
        Guide.yticks,
        "Component relative orientation")
    println("\r100%")
end


@step [run_rescaling_spread] function plot_rescaling_spread(p::Project)
    print("  0%")

    df = copy(p.step[run_rescaling_spread])
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

    ey = -180:15:180
    for i in 1:size(df, 1)
        for j in 2:length(ey)
            if df[i,:RelativeDir] < ey[j]
                df[i,:RelativeDir] = (ey[j-1] + ey[j]) / 2
                break
            end
        end
    end

    # Round scaling factors
    for scale in [1/2,1/√2,1,√2,2]
        df[df[:Scale].==scale,:Scale] = round(scale, 2)
    end
    print("\r 10%")

    function plot_col(col, tag, ymin, ymax, ylabel)
        dfm = by(df, [:RelativeDir, :Scale]) do d
            rho = col(d)
            m, s = mean(rho[!isnan(rho)]), sem(rho[!isnan(rho)])
            DataFrame(Mean=m, Min=m-s, Max=m+s)
        end

        h = plot(
            layer(dfm, x=:RelativeDir, ymin=:Min, ymax=:Max, color=:Scale, Geom.errorbar),
            layer(dfm, x=:RelativeDir, y=:Mean, color=:Scale, Geom.line, order=1),
            layer(dfm, x=:RelativeDir, y=:Mean, color=:Scale, Geom.point, order=2),
            Coord.cartesian(xmin=-180, xmax=180, ymin=ymin, ymax=ymax),
            Scale.color_discrete_manual(map(cmap, linspace(0, 0.9, 5))...),
            Scale.x_continuous(labels=x->@sprintf("%d°", x)),
            Guide.xticks(ticks=collect(-180:30:180)),
            Guide.xlabel("Angle relative to group direction"),
            Guide.ylabel(ylabel),
            Theme(; base_theme...))

        name = string("rescaling_", tag, ".pdf")
        draw(PDF(joinpath(plot_path, name), 6inch, 6inch/φ), h)
    end

    plot_col(d->d[:ExpectedCascadeSize],
        "cascade_size", 20, 70,
        "Expected cascade size")
    print("\r 40%")

    plot_col(d->d[:ExpectedCascadeSize] - d[:Detections],
        "spread", 0, 25,
        "Expected spread")
    print("\r 70%")

    plot_col(d->(d[:ExpectedCascadeSize] - d[:Detections]) ./ d[:Detections],
        "rel_spread", 0, 1.5,
        "Expected relative spread")
    println("\r100%")
end


@step [run_rescaling_spread_summary] function plot_rescaling_spread_summary(p::Project)
    print("  0%")

    df = p.step[run_rescaling_spread_summary]
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    print("\r 10%")

    h = plot(df, x=:Sensitivity, y=:Scale, color=:ExpectedCascadeSize,
        Geom.rectbin,
        Coord.cartesian(xmin=-1, xmax=1, ymin=-1, ymax=1),
        Scale.y_log2,
        Scale.color_continuous(colormap=cmap),
        Guide.xlabel("Sensitivity"),
        Guide.ylabel("Scaling factor"),
        Guide.colorkey("Expected\nCascade\nSize"),
        Theme(; base_theme...))
    name = string("rescaling_sensitivity_cascade_size.pdf")
    draw(PDF(joinpath(plot_path, name), 6inch, 6inch/φ), h)

    print("\r 40%")

    h = plot(df, x=:Sensitivity, y=:Scale, color=:ExpectedSpread,
        Geom.rectbin,
        Coord.cartesian(xmin=-1, xmax=1, ymin=-1, ymax=1),
        Scale.y_log2,
        Scale.color_continuous(colormap=cmap),
        Guide.xlabel("Sensitivity"),
        Guide.ylabel("Scaling factor"),
        Guide.colorkey("Expected\nSpread"),
        Theme(; base_theme...))
    name = string("rescaling_sensitivity_spread.pdf")
    draw(PDF(joinpath(plot_path, name), 6inch, 6inch/φ), h)

    print("\r 70%")

    h = plot(df, x=:Sensitivity, y=:Scale, color=:ExpectedRelativeSpread,
        Geom.rectbin,
        Coord.cartesian(xmin=-1, xmax=1, ymin=-1, ymax=1),
        Scale.y_log2,
        Scale.color_continuous(colormap=cmap),
        Guide.xlabel("Sensitivity"),
        Guide.ylabel("Scaling factor"),
        Guide.colorkey("Expected\nRelative\nSpread"),
        Theme(; base_theme...))
    name = string("rescaling_sensitivity_relative_spread.pdf")
    draw(PDF(joinpath(plot_path, name), 6inch, 6inch/φ), h)

    println("\r100%")
end



@step [run_rescaling_spread] function plot_rescaling_spread_sensitivity(p::Project)
    print("  0%")

    df = p.step[run_rescaling_spread]
    cmap = eval(p.conf[:colormap])
    base_theme = p.conf[:theme]

    # Limit to polarized
    df = df[df[:State].=="Polarized",:]

    print("\r  5%")

    # Limit to common disk
    rmax = minimum(df[:MaxRadius])
    df = df[df[:DistFromEdge] .<= rmax,:]

    print("\r 10%")

    # df[:Detections] = convert(DataVector{Float64}, df[:Detections])

    # Round scaling factors
    for scale in [1/2,1/√2,1,√2,2]
        df[df[:Scale].==scale,:Scale] = round(scale, 2)
    end
    print("\r 30%")

    dfm = by(df, [:Scale, :Sensitivity]) do d
        rho = d[:ExpectedCascadeSize]
        DataFrame(Mean=mean(rho[!isnan(rho)]))
    end

    print("\r 50%")

    h = plot(dfm, x=:Sensitivity, y=:Scale, color=:Mean,
        Geom.rectbin,
        Coord.cartesian(xmin=0, xmax=1, ymin=-1, ymax=1),
        Scale.y_log2,
        Guide.xlabel("Sensitivity"),
        Guide.ylabel("Scaling factor"),
        Theme(; base_theme...))

    print("\r 90%")

    name = string("rescaling_sensitivity_cascade_size.pdf")
    draw(PDF(joinpath(plot_path, name), 6inch, 6inch/φ), h)

    println("\r100%")
end


function runscale()
    lab = ["1÷2", "1÷√2", "1", "√2", "2"]
    val = [1/2, 1/√2, 1, √2, 2]
    println("> data/scale/fullscaled_$(lab[1]).h5")
    df = run_info("data/scale/fullscaled_$(lab[1]).h5", α = -0.2 / val[1])
    df[:Scale] = val[1]
    for i in 2:5
        println("> data/scale/fullscaled_$(lab[i]).h5")
        df2 = run_info("data/scale/fullscaled_$(lab[i]).h5", α = -0.2 / val[i])
        df2[:Scale] = val[i]
        append!(df, df2)
    end
    writetable("data/scale/df.csv", df)
    df
end

function runsize(mode)
    println("> size/data_5.0_$(mode).h5")
    _, _, dfb, dfe = run("size/data_5.0_$(mode).h5")
    dfb[:Size] = 5.0
    dfe[:Size] = 5.0
    for size in 10.0:5.0:25.0
        println("> size/data_$(size)_$(mode).h5")
        _, _, dfb2, dfe2 = run("size/data_$(size)_$(mode).h5", -0.2)
        dfb2[:Size] = size
        dfe2[:Size] = size
        append!(dfb, dfb2)
        append!(dfe, dfe2)
    end
    writetable("size/dfb_$(mode).csv", dfb)
    writetable("size/dfe_$(mode).csv", dfe)
    dfb, dfe
end


# Plotting
# --------

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

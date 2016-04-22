# Analyse effect of changing interindividual distances.

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

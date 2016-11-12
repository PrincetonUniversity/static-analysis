# Analyze effect of turbidity on schooling behavior.

function turbidity_project()
    Project("turbidity",
        α = -0.2,
        colormap = :plasma,
        theme = Dict(
            :minor_label_font => "CMU Serif",
            :minor_label_font_size => 10pt,
            :major_label_font => "CMU Serif",
            :major_label_font_size => 12pt,
            :key_title_font => "CMU Serif",
            :key_title_font_size => 12pt,
            :key_label_font => "CMU Serif",
            :key_label_font_size => 11pt,
            :plot_padding => 8pt,
        ),
    )
end

function turbidity_projects()
    p = Project[]
    for mode in ["control", "passive", "active"], λ in [1.0, 2.0, 4.0, 8.0], rep in 1:10
        push!(p, Project(format("turbidity_{1}_{2:.2f}_{3}", mode, λ, rep),
            mode = mode,
            λ = λ,
            replicate = rep,
            file = joinpath(data_path, format("turbidity_{1}_{2:.2f}_{3}.h5", mode, λ, rep)),
            α = -0.2,
            colormap = :plasma,
            theme = Dict(
                :minor_label_font => "CMU Serif",
                :minor_label_font_size => 10pt,
                :major_label_font => "CMU Serif",
                :major_label_font_size => 12pt,
                :key_title_font => "CMU Serif",
                :key_title_font_size => 12pt,
                :key_label_font => "CMU Serif",
                :key_label_font_size => 11pt,
                :plot_padding => 8pt,
            ),
        ))
    end
    p
end


@step function run_turbidity_info(p::Project)
    tps = turbidity_projects()
    for tp in tps
        println(format("\n# Turbidity: {} {:.2f}", tp.conf[:mode], tp.conf[:λ]))
        run_turbidity(tp)
        tp.step[run_turbidity][:Mode] = tp.conf[:mode]
        tp.step[run_turbidity][:AttenuationLength] = tp.conf[:λ]
    end
    vcat(DataFrame[tp.step[run_turbidity] for tp in tps])
end

ngroups(groups::Matrix{Int32}) = squeeze(maximum(groups, 1), 1)


function load_groups()
    df = DataFrame(
        [Int, ASCIIString, Float64, Int, Int32, Int, Int, Float64],
        [:Step, :Mode, :AttenuationLength, :Replicate, :GroupCount, :FissionRate, :FusionRate, :FissionFusionRate],
        0)
    for mode in ["control", "passive", "active"], λ in [1:.5:4;5:9], rep in 1:10
        path = joinpath(data_path, "turbidity", format("turbidity_{1}_{2:.1f}_{3}.h5", mode, λ, rep))
        groups = h5read(path, "groups")
        n, p = size(groups)
        count = ngroups(groups)
        dif = diff(count)
        fission = 10*vcat(0, dif.*(dif.>0))
        fusion = 10*vcat(0, -dif.*(dif.<0))
        fissionfusion = squeeze(mean([fission fusion], 2), 2)
        append!(df, DataFrame(Step=collect(1:p), Mode=mode,
            AttenuationLength=λ,
            Replicate=rep,
            GroupCount=count,
            FissionRate=fission,
            FusionRate=fusion,
            FissionFusionRate=fissionfusion))
    end
    df
end

function plot_groups(df::DataFrame)
    theme = Dict(
        :minor_label_font => "CMU Serif",
        :minor_label_font_size => 10pt,
        :major_label_font => "CMU Serif",
        :major_label_font_size => 12pt,
        :key_title_font => "CMU Serif",
        :key_title_font_size => 12pt,
        :key_label_font => "CMU Serif",
        :key_label_font_size => 11pt,
        :plot_padding => 8pt,
    )

    da = by(df, [:Mode,:AttenuationLength]) do d
        sub = d[:Step] .>= 1000
        m, s = mean(d[sub,:GroupCount]), std(d[sub,:GroupCount])
        DataFrame(MeanGroupCount=m, MinGroupCount=m-s, MaxGroupCount=m+s)
    end
    a = plot(da, x=:AttenuationLength, y=:MeanGroupCount, ymin=:MinGroupCount, ymax=:MaxGroupCount, color=:Mode, Geom.point, Geom.line, Geom.errorbar, Theme(; theme...), Guide.xlabel("Attenuation length (body lengths)"), Guide.ylabel("Average number of groups"))

    db = by(df, [:Mode,:AttenuationLength]) do d
        sub = d[:Step] .>= 1000
        m, s = mean(d[sub,:FissionRate]), sem(d[sub,:FissionRate])
        DataFrame(MeanFissionRate=m, MinFissionRate=m-s, MaxFissionRate=m+s)
    end
    b = plot(db, x=:AttenuationLength, y=:MeanFissionRate, ymin=:MinFissionRate, ymax=:MaxFissionRate, color=:Mode, Geom.point, Geom.line, Geom.errorbar, Theme(; theme...), Guide.xlabel("Attenuation length (body lengths)"), Guide.ylabel("Average fission rate (s<sup>-1</sup>)"))

    dc = by(df, [:Mode,:AttenuationLength]) do d
        sub = d[:Step] .>= 1000
        m, s = mean(d[sub,:FusionRate]), sem(d[sub,:FusionRate])
        DataFrame(MeanFusionRate=m, MinFusionRate=m-s, MaxFusionRate=m+s)
    end
    c = plot(dc, x=:AttenuationLength, y=:MeanFusionRate, ymin=:MinFusionRate, ymax=:MaxFusionRate, color=:Mode, Geom.point, Geom.line, Geom.errorbar, Theme(; theme...), Guide.xlabel("Attenuation length (body lengths)"), Guide.ylabel("Average fusion rate (s<sup>-1</sup>)"))

    dd = by(df, [:Mode,:AttenuationLength]) do d
        sub = d[:Step] .>= 1000
        m, s = mean(d[sub,:FissionFusionRate]), sem(d[sub,:FissionFusionRate])
        DataFrame(MeanFissionFusionRate=m, MinFissionFusionRate=m-s, MaxFissionFusionRate=m+s)
    end
    d = plot(dd, x=:AttenuationLength, y=:MeanFissionFusionRate, ymin=:MinFissionFusionRate, ymax=:MaxFissionFusionRate, color=:Mode, Geom.point, Geom.line, Geom.errorbar, Theme(; theme...), Guide.xlabel("Attenuation length (body lengths)"), Guide.ylabel("Average fission-fusion rate (s<sup><small>-1</small></sup>)"))

    a, b, c, d
end
# Config

data_path = joinpath(homedir(), "data", "static-analysis")
plot_path = joinpath(homedir(), "plots", "static-analysis")

base_project = Project("base",
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

function full_project()
    Project(base_project, "basic_full",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "fulldetectid.h5"),
        sensitivity = 1.0,
    )
end

function sub_project()
    Project(base_project, "basic_sub",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "subdetectid.h5"),
        sensitivity = 1.0,
    )
end

function full_blind_project()
    Project(base_project, "blind_full",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "fulldetectid.h5"),
    )
end

function random_projects()
    p = Project[]
    conf = Dict{AbstractString,Function}(
        "random_pos_vel_density" => makeRandomPosVelDensityDatasets,
        "random_pos_vel" => makeRandomPosVelDatasets,
        "random_pos_density" => makeRandomPosDensityDatasets,
        "random_pos" => makeRandomPosDatasets,
        "random_vel" => makeRandomVelDatasets,
        "control" => makeControlDatasets,
    )
    for (tag, f) in conf
        # f(joinpath(data_path, "fullprocessed.h5"), tag)
        push!(p, Project(base_project, tag,
            net_file = joinpath(data_path, "fullprocessed_$tag.h5"),
            detect_file = joinpath(data_path, "subdetectid_$tag.h5"),
            sensitivity = 1.0,
        ))
    end
    p
end

function random_project()
    Project(base_project, "basic_random",
        net_file = joinpath(data_path, "fullprocessed_random.h5"),
        detect_file = joinpath(data_path, "fulldetectid_random.h5"),
    )
end

function random_uniform_project()
    Project(base_project, "basic_random_uniform",
        net_file = joinpath(data_path, "fullprocessed_random_uniform.h5"),
        detect_file = joinpath(data_path, "subdetectid_random_uniform.h5"),
    )
end

function handpicked_hd_project()
    Project(base_project, "basic_handpicked_hd",
        net_file = joinpath(data_path, "handpickedprocessed.h5"),
        detect_file = joinpath(data_path, "handpickeddetectid_hd.h5"),
    )
end

function rescaling_project()
    Project(base_project, "rescaling_full",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "fulldetectid.h5"),
    )
end

function rescaled_projects()
    p = Project[]
    for s in [1/2, 1/√2, 1, √2, 2]
        push!(p, Project(base_project, format("rescaling_{:.2f}", s),
            scale = s,
            net_file = joinpath(data_path, format("fullprocessed_scale_{:.2f}.h5", s)),
            detect_file = joinpath(data_path, format("fulldetectid_scale_{:.2f}.h5", s)),
            α = -0.2 / s,
        ))
    end
    p
end

function rescaling_sub_project()
    Project(base_project, "rescaling_sub",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "subdetectid.h5"),
    )
end

function rescaled_sub_projects()
    p = Project[]
    for s in [1/2, 1/√2, 1, √2, 2], sensitivity in -1.0:0.2:1.0
        push!(p, Project(base_project, format("rescaling_sub_{:.2f}_{:.1f}", s, sensitivity),
            scale = s,
            sensitivity = sensitivity,
            net_file = joinpath(data_path, format("fullprocessed_scale_{:.2f}.h5", s)),
            detect_file = joinpath(data_path, format("subdetectid_scale_{:.2f}.h5", s)),
            α = -0.2 / s,
        ))
    end
    p
end
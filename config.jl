# Config

data_path = joinpath(homedir(), "data", "static-analysis")
plot_path = joinpath(homedir(), "plots", "static-analysis")

function full_project()
    Project("basic_full",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "fulldetectid.h5"),
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

function sub_project()
    Project("basic_sub",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "subdetectid.h5"),
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

function random_project()
    Project("basic_random",
        net_file = joinpath(data_path, "fullprocessed_random.h5"),
        detect_file = joinpath(data_path, "fulldetectid_random.h5"),
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

function rescaling_project()
    Project("rescaling_full",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "fulldetectid.h5"),
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

function rescaled_projects()
    p = Project[]
    for s in [1/2, 1/√2, 1, √2, 2]
        push!(p, Project(format("rescaling_{:.2f}", s),
            scale = s,
            net_file = joinpath(data_path, format("fullprocessed_scale_{:.2f}.h5", s)),
            detect_file = joinpath(data_path, format("fulldetectid_scale_{:.2f}.h5", s)),
            α = -0.2 / s,
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

function rescaling_sub_project()
    Project("rescaling_full",
        net_file = joinpath(data_path, "fullprocessed.h5"),
        detect_file = joinpath(data_path, "subdetectid.h5"),
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

function rescaled_sub_projects()
    p = Project[]
    for s in [1/2, 1/√2, 1, √2, 2]
        push!(p, Project(format("rescaling_sub_{:.2f}", s),
            scale = s,
            net_file = joinpath(data_path, format("fullprocessed_scale_{:.2f}.h5", s)),
            detect_file = joinpath(data_path, format("subdetectid_scale_{:.2f}.h5", s)),
            α = -0.2 / s,
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
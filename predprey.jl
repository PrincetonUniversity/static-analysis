# predprey model

using ODE
using DataFrames
using Gadfly
using Optim

odefun(a,b,c,d) = (t,x) -> [a*x[2]*(1-x[1])-b*x[1]; c*(1-x[1])*(1-x[2])-d*x[2]]

function run_risk()
    da = 0.5:0.5:9.5
    db = 1.0:3.0
    dc = 0.5:0.5:9.5
    dd = 1.0:2:9

    x0 = [0.0; 0.0]

    t = [0.0,50.0]

    # a = prey propensity to fragment
    # b = prey regrouping rate
    # c = pred propensity to coordinate
    # d = pred dispersion rate

    df = DataFrame([Float64, Float64, Float64, Float64, Float64, Float64, Float64],
        [:a, :b, :c, :d, :Risk, :Fragmentation, :Coordination], 0)

    for i in eachindex(da), j in eachindex(db), k in eachindex(dc), l in eachindex(dd)
        t, x = ode45(odefun(da[i], db[j], dc[k], dd[l]), x0, t, points=:specified)
        r = x[end][1]*(1-x[end][2])
        append!(df, DataFrame(a=da[i], b=db[j], c=dc[k], d=dd[l],
            Risk=r, Fragmentation=x[end][1], Coordination=x[end][2]))
    end

    df
end

function benefits(x)
    x0 = [0.0; 0.0]
    t = [0.0, 10.0]
    _, x = ode45(odefun(x[1], x[2], x[3], x[4]), x0, t, points=:specified)
    benef = x[end][1]*(1-x[end][2])
    coord = x[end][2]
    -min(benef, coord)
end

function run_opt(init)
    opt = minimize(Vector{Float64}(init))
    a, b, c, d = round(Optim.minimizer(opt), 2)
    x0 = [0.0; 0.0]
    t = [0.0, 10.0]
    _, x = ode45(odefun(a, b, c, d), x0, t, points=:specified)
    xinf, yinf = x[2]
    benef = xinf*(1-yinf)
    @printf("%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f\n", a, b, c, d, xinf, yinf, benef)
end

minimize() = minimize([1.0, 1.0, 1.0, 1.0])
function minimize(init)
    lower = [0.0, 0.0, 0.0, 0.0]
    upper = [Inf, Inf, Inf, Inf]
    optimize(DifferentiableFunction(benefits), init, lower, upper, Fminbox(), optimizer = GradientDescent)
end

function plot_risk(df)
    plot(df, y=:a, x=:c, ygroup=:b, xgroup=:d, color=:Risk,
        Scale.color_continuous(colormap=plasma),
        Guide.colorkey("x(1-y)"),
        Guide.ylabel("Prey propensity to fragment (a) <i><b>by</b></i> prey regrouping rate (b)"),
        Guide.xlabel("Predator propensity to coordinate (c) <i><b>by</b></i> predator dispersion rate (d)"),
        Geom.subplot_grid(Geom.rectbin),
        Theme(; Dict(
            :minor_label_font => "CMU Serif",
            :minor_label_font_size => 10pt,
            :major_label_font => "CMU Serif",
            :major_label_font_size => 12pt,
            :key_title_font => "CMU Serif",
            :key_title_font_size => 12pt,
            :key_label_font => "CMU Serif",
            :key_label_font_size => 11pt,
            :plot_padding => 8pt,
        )...))
end

nothing

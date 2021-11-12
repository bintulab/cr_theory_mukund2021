using Plots
pyplot(grid=false)
using DifferentialEquations
using ColorSchemes
using DataFrames
using CSV
using LaTeXStrings
using Statistics
using GLM
using Measures
using QuadGK
using Cairo, Fontconfig
import Gadfly
g = Gadfly
using Compose

fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

################################################################################
# Model:
# A <-> R <-> I
# ka: R -> A
# kr: A -> R
# ki: R -> I
# kp: I -> R
################################################################################

function three_state_model(dc, c, k, t)
    dc[1] = k['a'] * c[2] - k['s'] * c[1]
    dc[2] = k['s'] * c[1] - k['a'] * c[2] - k['i'] * c[2] + k['p'] * c[3]
    dc[3] = k['i'] * c[2] - k['p'] * c[3]
end

function isoutofdomain(c, k, t)
    if (c[1]<=1) && (c[2]<=1) && (c[3]<=1) && (c[1]>=0) && (c[2]>=0) && (c[3]>=0)
        return false
    else
        return true
    end
end

function run_recruitment(k_time_pairs, get_irrev=true)
    t = fill(Float64, 0)
    u = fill(Float64[], 0)

    c0 = [1.0; 0.0; 0.0]

    for pair in k_time_pairs
        k = pair[1]
        times = pair[2]
        prob = ODEProblem(three_state_model, c0, times, k)
        sol = solve(prob, isoutofdomain=isoutofdomain, saveat=0.01)
        t = [t; sol.t]
        u = [u; sol.u]
        c0 = sol.u[end]
    end
    as = [x[1] for x in u]
    rs = [x[2] for x in u]
    is = [x[3] for x in u]
    if get_irrev
        data = [as rs is]
    else
        data = [as rs]
    end

    return [t, data]
end

function plot_recruitment_2state(k_time_pairs, grey_pairs, blue_pairs, title, filename)
    sol = run_recruitment(k_time_pairs, false)
    p1 = Plots.plot()
    for pair in grey_pairs
        p1 = vspan!(p1, pair, label="", color=parse(Colorant, "#F5EBCE"))
    end
    for pair in blue_pairs
        p1 = vspan!(p1, pair, label="", color=parse(Colorant, "#DEE6F0"))
    end

    p1 = plot!(p1, 
               sol[1], sol[2],
               label=["Active" "Reversibly Silent"],
               color = [:green :grey],
               title=title,
               xlabel="Days",
               ylabel="Fraction of Cells",
               yticks = 0:0.5:1.0,
               xticks = [0, 5, 10, 15],
               linewidth=2,
               legendfont=12,
               legend=:topright,
               size=(500,250))

    savefig(p1, filename)
end

################################################################################
### FIGURE 1C
################################################################################

function plot_recruitment_2state_1c(k_time_pair_list, grey_pairs, blue_pairs, filename)
    sols = [run_recruitment(k, false) for k in k_time_pair_list]

    p1 = Plots.plot()
    for pair in grey_pairs
        p1 = vspan!(p1, pair, label="", color=parse(Colorant, "#F5EBCE"))
    end
    for pair in blue_pairs
        p1 = vspan!(p1, pair, label="", color=parse(Colorant, "#DEE6F0"))
    end

    colors = ["#009e73", "#0072b2"]
    labels = [L"k_s=10/\mathrm{d}, k_a=0.1/\mathrm{d}", L"k_s=0.1/\mathrm{d}, k_a=10/\mathrm{d}"]
    labels = [L"\mathrm{high} k_s\mathrm{, low} k_a", L"\mathrm{high} k_a\mathrm{, low} k_s"]
    labels = [L"high $k_s$, low $k_a$", L"low $k_s$, high $k_a$"]

    for i in 1:length(sols)
        sol = sols[i]
        p1 = plot!(p1, 
                   sol[1], sol[2][:,2],
                   xlabel="Days",
                   ylabel="Fraction Silent",
                   color = colors[i],
                   label = labels[i],
                   xlims = [0, 12],
                   ylims = [0, 1],
                   yticks = 0:0.5:1.0,
                   xticks = [0, 6, 12],
                   linewidth=2,
                   legendfont=12,
                   legend=:bottomright,
                   size=(500,250))
    end

    savefig(p1, filename)
end

k_mem_rec = Dict('a' => 0,
                 's' => 0.1, 
                 'i' => 0,
                 'p' => 0)
k_mem_rel = Dict('a' => 10,
                 's' => 0,
                 'i' => 0,
                 'p' => 0)
mem_time_pairs = [[k_mem_rec, (0.0, 6.0)],
                  [k_mem_rel, (6.0, 12.0)]]

k_sil_rec = Dict('a' => 0,
                 's' => 10,
                 'i' => 0,
                 'p' => 0)
k_sil_rel = Dict('a' => 0.1,
                 's' => 0,
                 'i' => 0,
                 'p' => 0)
sil_time_pairs = [[k_sil_rec, (0.0, 6.0)],
                  [k_sil_rel, (6.0, 12.0)]]

time_pairs = [sil_time_pairs, mem_time_pairs]

gp1 = [[0.0, 6.0]]
bp1 = [[6.0, 12.0]]
f1 = "fig_1c.pdf"

plot_recruitment_2state_1c(time_pairs, gp1, bp1, f1)

################################################################################
#### FIGURE 1D
################################################################################

c1 = parse(Colorant, "#1b9e77")
c2 = parse(Colorant, "#7570b3")
c3 = parse(Colorant, "#d95f02")

function cell_div_model(dc, c, k, t)
    ks = k['s']
    ka = k['a']
    kd = k['d']

    # number of on cells grows as gr, plus any reactivated off cells
    dc[1] = kd * c[1] + c[2] * ka - c[1] * ks
    # number of off cells grows as g, plus any silenced on cellsr
    dc[2] = kd * c[2] - c[2] * ka + c[1] * ks
end

function run_celldiv_model(ks, ka, kd, a0, tlim)
    c0 = [a0; 0.0] 

    println(c0)

    k = Dict('s' => ks,
             'a' => ka,
             'd' => kd)

    times = [0.0, tlim]
    prob = ODEProblem(cell_div_model, c0, times, k)
    sol = solve(prob, saveat=0.01)
    t = sol.t
    u = sol.u
    ons = [x[1] for x in u]
    offs = [x[2] for x in u]
    data = [offs ons]
    return [t, data]
end

function plot_celldiv(ks, ka, kd, a0, tlim)
    t, data = run_celldiv_model(ks, ka, kd, a0, tlim)
    p1 = Plots.plot()
    p1 = plot!(p1, 
               t, data,
               label = ["Reversibly Silent" "Active"],
               color = [c3 c2],
               xlabel = "Days",
               ylabel = "Number of Cells",
               ylims = [0, 700],
               yticks = [0, 350, 700],
               xlims = [0, 2],
               xticks = [0, 1, 2],
               linewidth = 2,
               legendfont = 12,
               size = (450, 300))
    return p1
end

p = plot_celldiv(1.0, 0.0, 1.0, 100, 2)
savefig(p, "fig_1d.pdf")

################################################################################
#### FIGURE 1E
################################################################################

signals = zeros(0)

function sinusoid(dc, c, k, t)
    dc[1] = k['a'] * c[2] - k['s'] * c[1]
    dc[2] = k['s'] * c[1] - k['a'] * c[2]
end
    

function isoutofdomain_sinus(c, k, t)
    for i in c
        if (i > 1) || (i < 0)
            return true
        end
    end
    return false
end

function run_sinus(k_time_pairs, get_irrev=true)
    t = fill(Float64, 0)
    u = fill(Float64[], 0)

    c0 = [0.75; 0.25]

    for pair in k_time_pairs
        k = pair[1]
        times = pair[2]
        prob = ODEProblem(sinusoid, c0, times, k)
        sol = solve(prob, isoutofdomain=isoutofdomain_sinus, saveat=0.001)
        t = [t; sol.t]
        u = [u; sol.u]
        c0 = sol.u[end]
    end
    as = [x[1] for x in u]
    rs = [x[2] for x in u]

    data = [as rs]
    return [t, data]
end

function get_plot_sinus(k_time_pairs, grey_pairs, blue_pairs, duration)
    sol = run_sinus(k_time_pairs)
    tlim = duration

    p1 = Plots.plot()

    crec = parse(Colorant, "#F5EBCE")
    crel = parse(Colorant, "#DEE6F0")

    for pair in grey_pairs
        ymin = 0.005
        ymax = 1.05
        p1 = plot!(p1, 
                   [pair[1], pair[1], pair[2], pair[2]],
                   [ymin, ymax, ymax, ymin],
                   color = crec,
                   linealpha = 0,
                   label = "",
                   seriestype = :shape)
        #  p1 = vspan!(p1, pair, label="", linecolor=crec, fillcolor=crec, linealpha=0, fillalpha=0)
    end
    for pair in blue_pairs
        ymin = 0.005
        ymax = 1.05
        p1 = plot!(p1, 
                   [pair[1], pair[1], pair[2], pair[2]],
                   [ymin, ymax, ymax, ymin],
                   color = crel,
                   linealpha = 0,
                   label = "",
                   seriestype = :shape)
        #  p1 = vspan!(p1, pair, label="", linecolor=crel, fillcolor=crel, linealpha=0)
    end

    p1 = plot!(p1, 
               sol[1], sol[2][:,2],
               order = 10,
               label="",
               legend=:bottomright,
               xlabel="Days",
               #  xticks= [0, 1, 2, 3, 5, 10],
               xticks = [0, tlim/2.0, tlim],
               xlims = (0, tlim),
               ylabel="Fraction Silent",
               yticks = [0, 0.5, 1],
               ylims = (0, 1.05),
               linewidth=2,
               color=:green,
               bottom_margin=0.25cm,
               right_margin=0.25cm,
               titlefontfamily = "Arial",
               size=(450,125))

    #  tf2_vals = sol[2][:, 1]
    #  tf2_vals_trunc = tf2_vals[Int(floor(0.5*length(tf2_vals))):length(tf2_vals)]

   
    return p1
end

function get_sinus_timepairs(ka, ks, period, duration)
    k_rec = Dict('a' => ka,
                 's' => ks)
    k_rel = Dict('a' => ka,
                 's' => 0)
    l = [[[k_rec, (i, i+period)], [k_rel, (i+period, i+2*period)]] for i in 0:2*period:duration]
    return [l[i][j] for i in 1:length(l) for j in 1:length(l[i])]
end

function get_greypairs(period, duration)
    return [[i, i + period] for i in 0:2*period:duration]
end

function get_bluepairs(period, duration)
    return [[i, i + period] for i in period:2*period:duration]
end

function get_plotpairs(period, ks)
    duration = 12
    ka = 1.0
    return [get_sinus_timepairs(ka, ks, period, duration),
            get_greypairs(period, duration),
            get_bluepairs(period, duration),
            duration]
end

palette = ColorSchemes.Spectral_4
palette = [colorant"#e66101", colorant"#fdb863", colorant"#b2abd2", colorant"#5e3c99"]
palette = ["grey25", "brown4", "blue", "darkviolet"]
palette = [colorant"skyblue", colorant"purple", colorant"green", colorant"darkorange"]
palette = [colorant"grey25", colorant"darkorange", colorant"purple", colorant"green"]
p1 = get_plot_sinus(get_plotpairs(3.0, 1.0)...)
p1 = plot!(p1, title=L"Long Period, Low $k_s$", titlefontcolor = palette[1], titleweight = "bold")

p2 = get_plot_sinus(get_plotpairs(0.25, 1.0)...)
p2 = plot!(p2, title=L"Short Period, Low $k_s$", titlefontcolor = palette[2])

p3 = get_plot_sinus(get_plotpairs(3.0, 3.0)...)
p3 = plot!(p3, title=L"Long Period, High $k_s$", titlefontcolor = palette[3])

p4 = get_plot_sinus(get_plotpairs(0.25, 3.0)...)
p4 = plot!(p4, title=L"Short Period, High $k_s$", titlefontcolor = palette[4])

lay = @layout [a b ; c d]
p = plot(p4, p3, p2, p1, layout=lay, size=(600, 450), linewidth=2)
savefig(p, "fig_1e.pdf")


################################################################################
#### FIGURE 1F
################################################################################

function bot_val(t, ks, ka)
    num = -exp(ka*t)*ka + exp(t*(-ka-ks))*ka + exp(t*(-ka-ks))*ks - exp(ka*t + t*(-ka - ks))*ks
    den = (exp(ka*t) - exp(t*(-ka - ks)))*(ka + ks)
    # den = (exp(t*(2*ka+ks)) - 1)*(ka + ks)
    return -num/den
end

function top_val(t, ks, ka)
    num = exp(ka*t)*ka - exp(t*(-ka-ks))*ka - ks + exp(ka*t)*ks
    den = (exp(ka*t) - exp(t*(-ka - ks)))*(ka + ks)
    return num/den
end

t_vals = zeros(0)
k_vals = zeros(0)
d_vals = zeros(0)

for t in 0:0.01:10
    for k in 0:0.1:10
        d = top_val(t/2.0, k, 1.0) - bot_val(t/2.0, k, 1.0)
        if t == 0
            d = 0.0
        end
        push!(t_vals, t)
        push!(k_vals, k)
        push!(d_vals, d)
    end
end

data = DataFrame(Dict("Period of Oscillation" => t_vals,
                      "Rate of Silencing (1/days)" => k_vals,
                      "Amplitude" => d_vals))

function get_pt_annotation(x, y, color="red")
    return g.Guide.annotation(compose(context(), 
                                      Compose.circle(x, y, 0.12), 
                                      Compose.stroke("white"),
                                      Compose.linewidth(0.5),
                                      fill(color), 
                                     ))
end
palette = [colorant"#e66101", colorant"#fdb863", colorant"#b2abd2", colorant"#5e3c99"]
palette = ["grey25", "brown4", "blue", "darkviolet"]
palette = [colorant"skyblue", colorant"purple", colorant"green", colorant"darkorange"]
palette = [colorant"grey25", colorant"darkorange", colorant"green", colorant"purple"]

p = g.plot(data, x="Period of Oscillation", y="Rate of Silencing (1/days)", color="Amplitude",
           g.Geom.rectbin, g.Guide.xticks(ticks=collect(0:1:5)), g.Guide.yticks(ticks=collect(0:1:5)),
           g.Coord.cartesian(xmin=0, xmax=5, ymin=0, ymax=5, aspect_ratio=1),
           g.Scale.ContinuousColorScale(p -> get(ColorSchemes.viridis, p)),
           get_pt_annotation(3.0, 1.0, palette[1]),
           get_pt_annotation(0.25, 1.0, palette[2]),
           get_pt_annotation(3.0, 3.0, palette[4]),
           get_pt_annotation(0.25, 3.0, palette[3]),
           g.Theme(major_label_font="ArialMT", minor_label_font="ArialMT",
                   key_title_font="ArialMT", key_label_font="ArialMT", key_label_color="black",
                   grid_color="black",
                   minor_label_color="black", major_label_color="black", key_title_color="black"))

img = PDF("fig_1f.pdf")
draw(img, p)


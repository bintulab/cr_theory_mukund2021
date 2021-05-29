using DifferentialEquations
using DataFrames
using CSV
using LaTeXStrings
using Statistics
using GLM
using StatsPlots
using Optim
using ForwardDiff
using ColorSchemes
using Gadfly
using Cairo
using Fontconfig
g = Gadfly

#
pyplot(grid=false)
#
fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

crec = parse(Colorant, "#F5EBCE")
crel = parse(Colorant, "#DEE6F0")

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

function run_recruitment(k_time_pairs, get_irrev=true, c0=[1.0; 0.0; 0.0])
    t = fill(Float64, 0)
    u = fill(Float64[], 0)

    for pair in k_time_pairs
        k = pair[1]
        times = pair[2]
        prob = ODEProblem(three_state_model, c0, times, k)
        sol = solve(prob, isoutofdomain=isoutofdomain, saveat=0.001)
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

################################################################################
### S2C-D
################################################################################

function get_pulsed_silencing_curve(s=200)
    pulse_pairs = [[k_recruit, (0.0, 1.0)],
                   [k_release, (1.0, 2.0)],
                   [k_recruit, (2.0, 3.0)],
                   [k_release, (3.0, 4.0)],
                   [k_recruit, (4.0, 5.0)],
                   [k_release, (5.0, 6.0)],
                   [k_recruit, (6.0, 7.0)],
                   [k_release, (7.0, 8.0)],
                   [k_recruit, (8.0, 9.0)],
                   [k_release, (9.0, 10.0)],
                   [k_release, (10.0, 15.0)]
                   ]
    gpp = [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0], [8.0, 9.0]]
    bpp = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [9.0, 15.0]]

    psol = run_recruitment(pulse_pairs)

    pp = plot()
    for pair in gpp 
        ymin = 0.005
        ymax = 1.05
        pp = plot!(pp, 
                   [pair[1], pair[1], pair[2], pair[2]],
                   [ymin, ymax, ymax, ymin],
                   color = crec,
                   linealpha = 0,
                   linewidth = 0,
                   label = "",
                   seriestype = :shape)
        #  pp = vspan!(pp, pair, label="", linecolor=crec, fillcolor=crec, linealpha=0, fillalpha=0)
    end
    for pair in bpp 
        ymin = 0.005
        ymax = 1.05
        pp = plot!(pp, 
                   [pair[1], pair[1], pair[2], pair[2]],
                   [ymin, ymax, ymax, ymin],
                   color = crel,
                   linealpha = 0,
                   linewidth = 0,
                   label = "",
                   seriestype = :shape)
        #  pp = vspan!(pp, pair, label="", linecolor=crel, fillcolor=crel, linealpha=0)
    end
    #  for pair in gpp
    #      pp = vspan!(pp, pair, label="", color=crec)
    #  end
    #  for pair in bpp
    #      pp = vspan!(pp, pair, label="", color=crel)
    #  end

    pp = plot!(pp,
               psol[1], psol[2],
               label=["A" "R" "I"],
               color = [:blue :orange :red],
               xlabel="Time (days)",
               ylabel="Fraction of Cells",
               yticks = 0:0.5:1.0,
               xticks = 0:5:15,
               #  xlims = (-Inf, 15),
               #  ylims=(-Inf, 1.0),
               linewidth=4,
               legend = :outerright,
               order = 0,
               size=(500, s))
    return pp
end

function get_pulsed_vs_continuous(s=200)
    pulse_pairs = [[k_recruit, (0.0, 1.0)],
                   [k_release, (1.0, 2.0)],
                   [k_recruit, (2.0, 3.0)],
                   [k_release, (3.0, 4.0)],
                   [k_recruit, (4.0, 5.0)],
                   [k_release, (5.0, 6.0)],
                   [k_recruit, (6.0, 7.0)],
                   [k_release, (7.0, 8.0)],
                   [k_recruit, (8.0, 9.0)],
                   [k_release, (9.0, 10.0)],
                   [k_release, (10.0, 15.0)]
                   ]
    gpp = [[0.0, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0], [8.0, 9.0]]

    cont_pairs = [[k_release, (0.0, 5.0)],
                  [k_recruit, (5.0, 10.0)],
                  [k_release, (10.0, 15.0)]]
    gpc = [[5.0, 10.0]]

    psol = run_recruitment(pulse_pairs)
    csol = run_recruitment(cont_pairs)

    pc = plot(
              psol[1], psol[2][:, 3],
              label="Pulsed",
              color=:purple,
              xlabel="Time (days)",
              ylabel="I(t)",
              ylims=(0.0, 1.0),
              yticks = 0:0.5:1.0,
              linewidth=3,
              legend = :topright,
              size=(500, 300))
    pc = plot!(pc,
               csol[1], csol[2][:, 3],
               label="Continuous",
               color=:deepskyblue,
               linewidth=3)


    for pair in gpp
        pc = vspan!(pc, pair, label="", color=:violet, alpha=0.2)
    end
    for pair in gpc
        pc = vspan!(pc, pair, label="", color=:cyan, alpha=0.2)
    end
    return pc
end

ksi = 1
ksa = 1

ks = 1

k_recruit = Dict('a' => ks/ksa,
                 's' => ks,
                 'i' => ks/ksi,
                 'p' => 0)

k_release = Dict('a' => ks/ksa,
                 's' => 0,
                 'i' => 0,
                 'p' => 0)

pp = get_pulsed_silencing_curve(500)
pc = get_pulsed_vs_continuous(500)

savefig(pp, "fig_s2c.pdf")
savefig(pc, "fig_s2d.pdf")

################################################################################
### S2E
################################################################################

function get_max_delta(ksa, ksi)
    mdelta = 0
    k_recruit = Dict('a' => 0.3,
                     's' => 0.3 * ksa,
                     'i' => 0.3 * ksa / ksi,
                     'p' => 0)
    k_release = Dict('a' => 0.3,
                     's' => 0,
                     'i' => 0,
                     'p' => 0)
    for t in collect(0.01:0.01:4.99)
        kt_pairs = [[k_recruit, (0.0, t)],
                    [k_release, (t, 5.0)]]
        sol_on = run_recruitment(kt_pairs)
        sol_off = run_recruitment(kt_pairs, true, [0.0; 1.0; 0.0])
        vals_on = sol_on[2][:,3]
        vals_off = sol_off[2][:,3]
        deltas = vals_off - vals_on
        md = maximum(deltas)
        mdelta = maximum([md, mdelta])
    end
    return mdelta
end

function analytical_max_delta(ka, ks, ki)
    function irrev_delta(t, ka, ks, ki)
        term1 = -0.5 * (ka + ki + ks + sqrt(-4 * ka * ki + (ka + ki + ks)^2)) * t
        term2 = (-1 + exp(sqrt(-4 * ka * ki + (ka + ki + ks)^2) * t)) * ki
        term3 = sqrt(-4 * ka * ki + (ka + ki + ks)^2)
        return (exp(term1) * term2)/term3
    end

    function tmax_delta(ka, ks, ki)
        term1 = -ka - ki - ks - sqrt(-4 * ki * ks + (ka + ki + ks)^2)
        term2 = -ka - ki - ks + sqrt(-4 * ki * ks + (ka + ki + ks)^2)
        term3 = sqrt(-4 * ki * ks + (ka + ki + ks)^2)
        return log(term1/term2)/term3
    end

    return irrev_delta(tmax_delta(ka, ks, ki), ka, ks, ki)
end

ksa_vals = zeros(0)
ksi_vals = zeros(0)
md_vals = zeros(0)

for ksa in 0.01:0.1:10.01
    for ksi in 0.01:0.1:10.01
        md = analytical_max_delta(1.0, 1.0 * ksa, 1.0 * ksa/ksi)
        push!(ksa_vals, ksa)
        push!(ksi_vals, ksi)
        push!(md_vals, md)
    end
end

data = DataFrame(Dict("Ratio of Reversible Silencing to Activation" => ksa_vals, 
                      "Ratio of Reversible to Irreversible Silencing" => ksi_vals, 
                      "Gap" => md_vals))

println("plotting")
p = g.plot(data, 
         x="Ratio of Reversible to Irreversible Silencing",
         y="Ratio of Reversible Silencing to Activation", 
         color="Gap", g.Geom.rectbin,
         g.Guide.yticks(ticks=collect(0:10)), g.Guide.xticks(ticks=collect(0:10)),
         g.Guide.xlabel("ks/ka ratio"), g.Guide.ylabel("ks/ki ratio"),
         g.Coord.cartesian(xmin=0, xmax=10, ymin=0, ymax=10, aspect_ratio=1),
         g.Scale.ContinuousColorScale(p -> get(ColorSchemes.viridis, p)),
         g.Theme(major_label_font="ArialMT", minor_label_font="ArialMT",
               key_title_font="ArialMT", key_label_font="ArialMT", key_label_color="black",
               grid_color="black", 
               minor_label_color="black", major_label_color="black", key_title_color="black"))

img = PDF("fig_s2e.pdf")
draw(img, p)

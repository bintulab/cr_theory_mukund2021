using DifferentialEquations, Random, Distributions
using Plots, ColorSchemes
using DataFrames, CSV, LaTeXStrings
using Statistics, GLM, Measures
using LightGraphs, SimpleWeightedGraphs#, GraphPlot
using Compose, Cairo, Fontconfig
using LsqFit, Dates
import Gadfly
g = Gadfly # so we can use gadfly too! 

pyplot(grid=false)
#  plotly()

epsilon = 1e-8

fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

################################################################################
#### FIGURE 2B
################################################################################

c_starts = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

function one_cr_model(dc, c, k, t)
    ka = k['a']
    ks = k['s']
    ki = k['i']

    dc[1] = ka * c[2] - ks * c[1] 
    dc[2] = ks * c[1] - ka * c[2] - ki * c[2] * c[1]
    dc[3] = ki * c[2] * c[1]

    @assert abs(sum(dc))<=epsilon
end

function run_one_cr_recruitment(c0_val)
    t = fill(Float64, 0)
    u = fill(Float64[], 0)

    k_sil = Dict('a' => 1.0,
                 's' => 1.0,
                 'i' => 0,
                 'p' => 0)
    kt_pairs = [[k_sil, (0.0, 2.0)]]

    c0 = [c0_val; 1.0-c0_val; 0.0]

    for pair in kt_pairs
        k = pair[1]
        times = pair[2]
        prob = ODEProblem(one_cr_model, c0, times, k)
        sol = solve(prob, saveat=0.001)
        t = [t; sol.t]
        u = [u; sol.u]
        c0 = sol.u[end]
    end
    a1s = [x[1] for x in u]

    return [t, a1s]
end

function plot_c0_vals(c0_vals)
    cscheme = get(ColorSchemes.viridis, range(0.0, stop=0.9, length=length(c0_vals)))
    p = plot()
    for i in 1:length(c0_vals)
        c0 = c0_vals[i]
        t, d = run_one_cr_recruitment(c0)
        p = plot!(p, t, d, 
                 label=string(c0), color=cscheme[i], 
                 xlabel="Days",
                 ylabel="Fraction Active",
                 yticks = 0:0.5:1.0,
                 ylims=(-0.1, 1.1),
                 xlims=(-0.1, 2.1),
                 xticks=0:1.0:2.0,
                 legend=:outerright,
                 legendtitle="A(0)",
                 #  legendtitlefontsize=14,
                 linewidth=2,
                 size=(600, 300))
    end

    return p
end

savefig(plot_c0_vals(c_starts), "fig_2b.pdf")

################################################################################
#### FIGURE 2C
################################################################################
function steady_state(r)
    ka = 1.0
    ks = r
    #  return (-ka + sqrt(ka^2 + 4*ka*ks))/(2*ks)
    return ka/(ka + ks)
end

xs = collect(0:0.001:10)
ys = [steady_state(x) for x in xs]

fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

p = plot(xs, ys,
         xlabel=L"$k_s/k_a$",
         ylabel="Steady-State 
Fraction Active",
         label="", color=:green,
         yticks = [0.0, 1.0],
         ylims=(-0.1, 1.1),
         xlims=(-0.1, 10),
         xticks=[0, 10],
         linewidth=2,
         size=(300, 261))
savefig(p, "fig_2c.pdf")

################################################################################
#### FIGURE 2E
################################################################################

function ksval(p)
    exp = 3.0
    kd = 0.33
    return (p^exp)/(p^exp + (kd)^exp) 
end

pvals = 0:0.01:100
kvals = [ksval(p) for p in pvals]

p = plot(pvals, kvals,
         linewidth=2,
         xlims=[0, 1.5],
         xticks=[0, 1],
         xformatter = x -> [L"$0$", L"$[CR]_{max}$"][Int(x+1)],
         ylims=[0, 1],
         yticks=[0, 1],
         yformatter = y -> [L"$0$", L"$k_{s\ max}$"][Int(y+1)],
         xlabel="[CR]",
         ylabel=L"$k_s($CR$)$",
         label=false,
         size=(300, 261))
savefig(p, "fig_2e.pdf")

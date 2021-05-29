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
using Random, Distributions
using Plots, ColorSchemes
using GLM, Measures
using LightGraphs, SimpleWeightedGraphs#, GraphPlot
using Compose, Cairo, Fontconfig
using LsqFit, Dates
import Gadfly
g = Gadfly 
#
pyplot(grid=false)
#
fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

crec = parse(Colorant, "#F5EBCE")
crel = parse(Colorant, "#DEE6F0")

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

function plot_recruitment(k_time_pairs, grey_pairs, blue_pairs, filename)
    sol = run_recruitment(k_time_pairs)

    p1 = plot()
    for pair in grey_pairs
        p1 = vspan!(p1, pair, label="", color=crec)
    end
    for pair in blue_pairs
        p1 = vspan!(p1, pair, label="", color=crel)
    end
    p1 = plot!(p1,
               sol[1], sol[2],
               label=["Active (A)" "Reversibly Silent (R)" "Irreversibly Silent (I)"],
               color = [:blue :orange :red],
               xlabel="Days",
               ylabel="Fraction of Cells",
               yticks = 0:0.2:1.0,
               linewidth=4,
               legendfontsize=13,
               size=(600,300))
    savefig(p1, filename)
end

################################################################################
### FIGURE 3A
################################################################################

k_recruit = Dict('a' => 0.1,
                 's' => 1,
                 'i' => 0.2,
                 'p' => 0)

k_release = Dict('a' => 0.1,
                 's' => 0,
                 'i' => 0,
                 'p' => 0)
kt_pairs = [[k_recruit, (0.0, 5.0)],
            [k_release, (5.0, 15.0)]]
gp = [[0.0, 5.0]]
bp = [[5.0, 15.0]]
fname = "fig_3a.pdf"

plot_recruitment(kt_pairs, gp, bp, fname)

################################################################################
### FIGURE 3B-C
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

function plot_pulsed_comparison(s=200)
    l = @layout [a ; b ; c]

    pp = get_pulsed_silencing_curve(s)
    pc = get_pulsed_vs_continuous(s)

    savefig(pp, "fig_3b.pdf")
    savefig(pc, "fig_3c.pdf")
    # p = plot(pp, pc, ps, layout=l)
    # f2 = "fig_2c.pdf"
    # savefig(p, f2)
end

plot_pulsed_comparison()


################################################################################
### FIGURE 3D
################################################################################


function three_state_model(dc, c, k, t)
    dc[1] = k['a'] * c[2] - k['s'] * c[1]
    dc[2] = k['s'] * c[1] - k['a'] * c[2] - k['i'] * c[2] + k['p'] * c[3]
    dc[3] = k['i'] * c[2] - k['p'] * c[3]

    @assert abs(dc[1]+dc[2]+dc[3])<=epsilon
end

function isoutofdomain(c, k, t)
    for i in 1:3:length(c)
        if abs(1-(c[i] + c[i+1] + c[i+2])) >= epsilon
            println(c)
            println(k)
            println(t)
            return false
        end
    end
    return true
end

function run_recruitment(k_time_pairs, c0=Nothing)
    t = fill(Float64, 0)
    u = fill(Float64[], 0)

    if (c0 == Nothing)
        c0 = [1.0; 0.0; 0.0]
    end

    for pair in k_time_pairs
        k = pair[1]
        times = pair[2]
        prob = ODEProblem(three_state_model, c0, times, k)
        sol = solve(prob, #AutoTsit5(Rosenbrock23());
                    # isoutofdomain=isoutofdomain,
                    saveat=0.01,
                    d_discontinuities=[times[1], times[2]])
        t = [t; sol.t]
        u = [u; sol.u]
        c0 = sol.u[end]
    end
    as = [x[1] for x in u]
    rs = [x[2] for x in u]
    is = [x[3] for x in u]
    data = [as rs is]

    return [t, data]
end

k_strong_sil_recruit = Dict(
    'a' => 0.31,
    's' => 11,
    'i' => 0.3,
    'p' => 0)

k_strong_sil_release = Dict(
    'a' => 0.31,
    's' => 0.1,
    'i' => 0,
    'p' => 0)

k_weaka_recruit = Dict(
    'a' => 0.69,
    's' => 0.1,
    'i' => 0,
    'p' => 0.1)

k_nomema_recruit = Dict('a' => 1.0,
                        's' => 0.1,
                        'i' => 0,
                        'p' => 0)

k_nomema_release = Dict('a' => 0,
                        's' => 0.1,
                        'i' => 0,
                        'p' => 0)

function plot_act(kp, kr, fname) # ka ks
    k_rec = Dict('a' => kp, 's' => kr, 'i' => 0, 'p' => 0)
    k_rel = Dict('a' => 0,  's' => kr, 'i' => 0, 'p' => 0)
    pair1 = [k_rec, (0, 5.0)]
    pair2 = [k_rel, (5.0, 15.0)]
    c0 = [0.0; 1.0; 0.0]
    sol = run_recruitment([pair1, pair2], c0)
    seq = plot()
    seq = vspan!(seq, [0, 5], label="", color=parse(Colorant, "#F5EBCE"))
    seq = vspan!(seq, [5, 15], label="", color=parse(Colorant, "#DEE6F0"))
    seq = plot!(seq, 
                sol[1], sol[2][:, 2],
                color = :grey,
                label = "",
                xlabel="Days",
                xlims=[0, 15],
                xticks = 0:5:15,
                ylabel="Fraction Silent",
                yticks = 0:0.5:1.0,
                ylims = [0, 1.0],
                borderpad = 0.0,
                legendfontsize=13,
                linewidth=2)

    #  seq.o[:legend](borderaxespad=-1.0)
    seq = plot!(seq, 
                size=(600, 300),
                legend=true)
    savefig(seq, fname)
    return seq
end

function plot_sil(ka, ks, ki, fname)
    k_rec = Dict('a' => ka, 's' => ks, 'i' => ki, 'p' => 0)
    k_rel = Dict('a' => ka, 's' => 0,  'i' => 0,  'p' => 0)
    pair1 = [k_rec, (0, 5.0)]
    pair2 = [k_rel, (5.0, 15.0)]
    c0 = [1.0; 0.0; 0.0]
    sol = run_recruitment([pair1, pair2], c0)
    seq = plot()
    seq = vspan!(seq, [0, 5], label="", color=parse(Colorant, "#F5EBCE"))
    seq = vspan!(seq, [5, 15], label="", color=parse(Colorant, "#DEE6F0"))
    seq = plot!(seq, 
                sol[1], [1.0-a for a in sol[2][:, 1]],
                color = :grey,
                label = "",
                xlabel="Days",
                xlims=[0, 15],
                xticks = 0:5:15,
                ylabel="Fraction Silent",
                yticks = 0:0.5:1.0,
                ylims=[0, 1.0],
                borderpad = 0.0,
                legendfontsize=13,
                linewidth=2)

    #  seq.o[:legend](borderaxespad=-1.0)
    seq = plot!(seq, 
                size=(600, 300),
                legend=true)
    savefig(seq, fname)
    return seq
end

plot_act(1.0, 1.0, "fig_3d_1.pdf")
plot_sil(1.1, 15.0, 0.5, "fig_3d_2.pdf")


################################################################################
### FIGURE 3E
################################################################################

function unpack_values(ari)
    return [ari['a'], ari['r'], ari['i']]
end

function unpack_parms(k)
    return [k['a'], k['s'], k['i'], k['p']]
end

function hill(k, a, exp=5, kd=0.5)
    if exp==Nothing
        exp = 5
    end
    if kd==Nothing
        kd = 0.5
    end
    return k * a #(a^exp/(a^exp + kd^exp))
end

function get_contrib(src_vals, dst_vals, kvals, he=Nothing, hk=Nothing, is_input=false, input_fn_val=Nothing)
    sa, sr, si     = src_vals
    da, dr, di     = dst_vals 
    ka, ks, ki, kp = unpack_parms(kvals) 
    
    if (sa == 0)
        return [0, 0, 0]
    end

    if (is_input)&&(input_fn_val != Nothing)
        sa = input_fn_val
    end

    dda =  dr*ka*sa - da*ks*sa
    ddr = -dr*ka*sa + da*ks*sa - dr*ki*sa + di*kp*sa
    ddi =                     dr*ki*sa - di*kp*sa

    # dda =  dr*ka - da*hill(ks, sa, he, hk)
    # ddr = -dr*ka + da*hill(ks, sa, he, hk) - dr*hill(ki, sa, he, hk) + di*hill(kp, sa, he, hk)
    # ddi =                                    dr*hill(ki, sa, he, hk) - di*hill(kp, sa, he, hk)
    
    @assert dda + ddr + ddi <= epsilon
    
    return [dda, ddr, ddi]
end

function get_bkgd_contrib(src_vals, kvals)
    sa, sr, si     = src_vals
    ka, ks, ki, kp = unpack_parms(kvals) 

    dda =  sr*ka - sa*ks
    ddr = -sr*ka + sa*ks - sr*ki + si*kp
    ddi =                  sr*ki - si*kp
    
    @assert dda + ddr + ddi <= epsilon
   
    return [dda, ddr, ddi]
end


function get_k_bkgd_func(a=0.0, s=0.0, i=0.0, p=0.0)
    return function()
        k_bkgd = Dict('a' => a,
                      's' => s,
                      'i' => i,
                      'p' => p)
        return k_bkgd
    end
end

function get_kval_func(k_sil_delta, k_act_delta)
    f = function get_k_vals(src, dst, weight, t, input_id, t_spans, kvspec)
        k_bkgd = Dict('a' => 0.0,
                      's' => 0.0,
                      'i' => 0.0,
                      'p' => 0.0)

        function add_silencing(k)
            k['s'] += k_sil_delta['s']
            k['i'] += k_sil_delta['i']
            return k
        end

        function add_activation(k)
            k['a'] += k_act_delta['a']
            k['p'] += k_act_delta['p']
            return k
        end

        function get_k(weight)
            if (weight==1)
                return add_activation(k_bkgd)
            elseif (weight==-1)
                return add_silencing(k_bkgd)
            end
        end

        # if src is the input node, and t is not within the timezones, k values are 0 
        # else return based on wghts 
        if (src == input_id) 
            input_on = false
            for span in t_spans
                t_start = span[1]
                t_end = span[2]
                if (t_start <= t) && (t <= t_end)
                    input_on = true
                end
            end
            if (! input_on) 
                return  k_bkgd
            end
        end

        if kvspec != Nothing
            kvs = get(kvspec, src, Nothing)
            if kvs != Nothing
                return kvs
            end
        end
        return get_k(weight)
    end
    return f
end

function digraph_model(du, u, p, t)
    digraph    = p["digraph"]
    input_id   = p["input_id"]
    t_spans    = p["t_spans"]
    hill_exp   = p["hill_exp"]
    hill_kd    = p["hill_kd"]
    get_k_vals = p["get_k_vals"]
    get_k_bkgd = p["get_k_bkgd"]
    kvspec = get(p, "k_spec", Nothing)
    input_fn = get(p, "input_fn", Nothing)

    if (input_fn != Nothing)
        input_fn_val = input_fn(t)
    else
        input_fn_val = Nothing
    end

    if get_k_vals == Nothing
        get_k_vals = get_kval_func(10.0, 1.0)
    end
    if get_k_bkgd == Nothing
        get_k_bkgd = get_k_bkgd_func()
    end

    vs = collect(vertices(digraph))
    es = collect(edges(digraph))
    # vertex ==> a, r, i
    vdict = Dict(v => Dict("da" => 0.0,
                           "dr" => 0.0,
                           "di" => 0.0) for v in vs)

    # compute deltas
    for v in vs 
        if (v != input_id)
            s_vals = u[3 * (v - 1) + 1 : (3 * (v - 1)) + 3] # compute background activation
            d_vals = s_vals
            k_vals = get_k_bkgd()
            deltas = get_bkgd_contrib(s_vals, k_vals)
            vdict[v]["da"] += deltas[1]
            vdict[v]["dr"] += deltas[2]
            vdict[v]["di"] += deltas[3]
            for e in [x for x in es if x.dst==v] # compute recruitment-driven activation/silencing
                src = e.src
                dst = e.dst

                @assert dst==v

                is_input = (src==input_id)

                weight = e.weight
                src_vals = u[(3 * (src - 1)) + 1 : (3 * (src - 1)) + 3]
                dst_vals = u[(3 * (dst - 1)) + 1 : (3 * (dst - 1)) + 3]
                k_vals = get_k_vals(src, dst, weight, t, input_id, t_spans, kvspec)
                deltas = get_contrib(src_vals, dst_vals, k_vals, hill_exp, hill_kd, is_input, input_fn_val)

                # if (deltas[1] > 0)&&(dst==length(vs))
                #     println(src)
                #     println(dst)
                #     println(weight)
                #     println(src_vals)
                #     println(dst_vals)
                #     println(k_vals)
                #     println(deltas)
                #     println()
                # end
                #
                vdict[v]["da"] += deltas[1]
                vdict[v]["dr"] += deltas[2]
                vdict[v]["di"] += deltas[3]
            end
        end
    end

    # return deltas inplace
    for idx in 1:length(vs)
        v = vs[idx]
        du[3 * (v-1) + 1] = vdict[v]["da"]
        du[3 * (v-1) + 2] = vdict[v]["dr"]
        du[3 * (v-1) + 3] = vdict[v]["di"]
    end
end

function run_digraph(parms, tlim, tstep=Nothing)
    digraph = parms["digraph"]
    num_vertices = length(collect(vertices(digraph)))
    
    if (parms["u0"] == Nothing)
        u0 = reduce(append!, [[1.0, 0.0, 0.0] for x in 1:num_vertices])
    else
        u0 = parms["u0"]
    end

    @assert length(u0) == 3*num_vertices

    tspan = (0.0, tlim)

    prob = ODEProblem(digraph_model, u0, tspan, parms)
    
    if (tstep==Nothing)
        sol = solve(prob, d_discontinuities=collect(Iterators.flatten(parms["t_spans"])))#, saveat=0.001)
    else
        sol = solve(prob, d_discontinuities=collect(Iterators.flatten(parms["t_spans"])), saveat=tstep)
    end

    t = sol.t
    u = sol.u

    avals = [[x[(3*(i-1))+1] for x in u] for i in 1:num_vertices]
    data = hcat(avals...)

    return [t, data]
end

function get_digraph_data(iff_parms, t_span)
    sol = run_digraph(iff_parms, t_span)
    labels = [string("CR ", i) for i in 1:length(collect(vertices(iff_parms["digraph"])))]
    labels = hcat(labels...)
    pdata = sol[2]
    return [sol[1], pdata, labels]
end

function plot_digraph(iff_parms, t_span, plot_input=true, plot_outputs_only=false, plot_last_only=false)
    sol = run_digraph(iff_parms, t_span)
    labels = [string("CR ", i) for i in 1:length(collect(vertices(iff_parms["digraph"])))]
    labels = hcat(labels...)
    if plot_input
        pdata = sol[2]
    elseif plot_outputs_only
        pdata = sol[2][:, end-1:end]
    elseif plot_last_only
        pdata = sol[2][:, end]
    else
        pdata = sol[2][:, 2:end]
    end
    p1 = plot(sol[1], pdata,
              label=labels,
              palette=:seaborn_colorblind,
              xlabel="Days",
              ylabel="Fraction Active",
              ylims=(-0.1, 1.1),
              yticks=0:0.5:1.0,
              xticks=0:t_span/2.0:t_span,
              linewidth=2,
              legend=:outerright,
              size=(700, 300))

    for span in iff_parms["t_spans"]
        p1 = vspan!(p1, span, label="", color=:gray, alpha=0.2)
    end
    return p1
end

function construct_parms(digraph, iid, spans, 
                          u0=Nothing, hill_exp=Nothing, hill_kd=Nothing,
                          get_k_bkgd=Nothing, get_k_vals=Nothing)
    parms = Dict("digraph"  => digraph,
                  "input_id" => iid,
                  "t_spans"  => spans,
                  "u0" => u0,
                  "hill_exp" => hill_exp,
                  "hill_kd" => hill_kd,
                  "get_k_vals" => get_k_vals,
                  "get_k_bkgd" => get_k_bkgd)
    return parms
end

colors = get(ColorSchemes.viridis, reverse([0.3, 0.6, 0.9]))

sources = [1, 1, 2]
destinations = [2, 3, 3]
wghts = [1, 1, -1]
digraph = SimpleWeightedDiGraph(sources, destinations, wghts)
iid = 1
spans = [[20, 40], [60, 80]]
u0 = [1.0, 0.0, 0.0, 
      0.0, 1.0, 0.0,
      0.0, 1.0, 0.0]
get_k_bkgd = get_k_bkgd_func(0, 1.0)

function ifn(t)
    if t <= 20
        return 0
    elseif t <= 40
        return 1
    elseif t<= 60 
        return 0
    elseif t<= 80
        return 1
    else
        return 0
    end
end

iff_parms = construct_parms(digraph, iid, spans, u0, Nothing, Nothing, get_k_bkgd, Nothing)
push!(iff_parms, "k_spec" => [Dict('a' => 1.0, 's' => 0.0,  'i' => 0, 'p' => 0), 
                               Dict('a' => 0.1, 's' => 15.0, 'i' => 0.5, 'p' => 0),
                               Dict('a' => 2.0, 's' => 0.0,  'i' => 0, 'p' => 0)])
# push!(iff_parms, "input_fn" => ifn)

colors  = ["turquoise4" "chocolate2" "rebeccapurple"]
colors = ColorSchemes.Dark2_3
t, data, ls = get_digraph_data(iff_parms, 100.0)
t_span = 100

p1 = plot(ifn, 0, 100,
          label="", 
          xlabel="", 
          ylabel="Input",
          guidefontsize = 18,
          ylims=(-0.1, 1.1),
          yticks=0:0.5:1.0,
          xticks=0:t_span/2.0:t_span,
          color=colors[1],
          linewidth=4)

p2 = plot(t, data[:, end-1], label = "", 
          xlabel="",
          ylabel="Neg. Arm",
          guidefontsize = 18,
          ylims=(-0.1, 1.1),
          yticks=0:0.5:1.0,
          xticks=0:t_span/2.0:t_span,
          linewidth=4,
          legend=:topright,
          color=colors[2],
          size=(600, 400))

p3 = plot(t, data[:, end], label = "", 
          xlabel="Days",
          ylabel="Output",
          guidefontsize = 18,
          ylims=(-0.02, 0.2),
          yticks=0:0.1:0.2,
          xticks=0:t_span/2.0:t_span,
          linewidth=4,
          legend=:topright,
          color=colors[3],
          size=(600, 500))

lay = @layout [a; b; c]
p = plot(p1, p2, p3, layout=lay, link=:all)
savefig(p, "fig_3e.pdf")

################################################################################
### FIGURE 3 + S2 Databases
################################################################################

function get_adapt_response(ka, ks, ki)
    sources = [1, 1, 2]
    destinations = [2, 3, 3]
    wghts = [1, 1, -1]
    digraph = SimpleWeightedDiGraph(sources, destinations, wghts)
    iid = 1
    spans = [[20, 40], [60, 80]]
    u0 = [1.0, 0.0, 0.0, 
          0.0, 1.0, 0.0,
          0.0, 1.0, 0.0]
    get_k_bkgd = get_k_bkgd_func(0, 1.0)
    iff_parms = construct_parms(digraph, iid, spans, u0, Nothing, Nothing, get_k_bkgd, Nothing)
    push!(iff_parms, "k_spec" => [Dict('a' => ka, 's' => 0.0,  'i' => 0, 'p' => 0), 
                                  Dict('a' => 0.1, 's' => ks, 'i' => ki, 'p' => 0),
                                  Dict('a' => 0.0, 's' => 0.0,  'i' => 0, 'p' => 0)])

    t, data, ls = get_digraph_data(iff_parms, 100.0)
    output = data[:, end]

    peak = maximum(output)
    pkind = findfirst(o->o==peak, output)

    sind = findfirst(time->time>=30, t)
    set = output[sind]

    @assert t[sind] >= t[pkind]
    @assert set >= 0
    @assert peak > 0

    k = -log(set/peak)/(t[sind] - t[pkind])

    pk2ind = findfirst(time->time>=50, t)
    pk2 = maximum(output[pk2ind:end])

    return [peak, log(2)/k, pk2]
end

function get_response_dbase_3f(small = false)
    karange = 0.01:0.1:20.01
    ksrange = 0.01:0.1:20.01
    kirange = [0.5]

    if small
        karange = 1:2
        ksrange = 1:10
        kirange = 0.25:0.25:0.5
    end

    dblen = length(karange) * length(ksrange) * length(kirange)
    println("Should take approx. " * string(dblen * 0.08) * " seconds")

    kas = zeros(dblen)
    kss = zeros(dblen)
    kis = zeros(dblen)
    pks = zeros(dblen)
    pk2s = zeros(dblen)
    ads = zeros(dblen)
    global idx = 1

    for ks in ksrange
        println("Working on ks ", ks)
        for ka in karange
            for ki in kirange
                pk, ad, pk2 = get_adapt_response(ka, ks, ki)
                kas[idx] = ka
                kss[idx] = ks
                kis[idx] = ki
                pks[idx] = pk
                pk2s[idx] = pk2
                ads[idx] = ad
                global idx = idx + 1
            end
        end
    end

    df = DataFrame(Dict("ka" => kas,
                        "ks" => kss,
                        "ki" => kis,
                        "ad" => ads,
                        "pk" => pks,
                        "pk2" => pk2s))
    if (! small) 
        CSV.write("../loop_simulations/3f_response_peaks.csv", df)
    end
end

function get_response_dbase(small = false)
    #  karange = sort([collect(0.01:0.5:20.01); collect(0.01:0.05:0.5); 0.5])
    #  ksrange = sort([collect(0.01:0.5:20.01); collect(0.01:0.05:0.5); 0.5])
    #  kirange = sort([collect(0.01:0.2:2.01); collect(0.01:0.02:0.21); 0.5])
    #
    #  # spanning in log space
    #  karange = sort([exp10.(range(log10(0.01), log10(20.01), length=20)); 0.5])
    #  ksrange = sort([exp10.(range(log10(0.01), log10(20.01), length=20)); 0.5])
    #
    karange = 0.01:1.0:20.01
    ksrange = 0.01:1.0:20.01

    karange = sort([collect(karange); 0.05; 0.10; 0.25; 0.5; 1.0])
    ksrange = sort([collect(ksrange); 0.05; 0.10; 0.25; 0.5; 1.0])
    kirange = sort([exp10.(range(log10(0.01), log10(2.01), length=length(karange))); 0.5])

    if small
        karange = 1:2
        ksrange = 1:10
        kirange = 0.25:0.25:0.5
    end

    dblen = length(karange) * length(ksrange) * length(kirange)
    println("Should take approx. " * string(dblen * 0.08) * " seconds")

    kas = zeros(dblen)
    kss = zeros(dblen)
    kis = zeros(dblen)
    pks = zeros(dblen)
    pk2s = zeros(dblen)
    ads = zeros(dblen)
    global idx = 1

    for ks in ksrange
        println("Working on ks ", ks)
        for ka in karange
            for ki in kirange
                pk, ad, pk2 = get_adapt_response(ka, ks, ki)
                kas[idx] = ka
                kss[idx] = ks
                kis[idx] = ki
                pks[idx] = pk
                pk2s[idx] = pk2
                ads[idx] = ad
                global idx = idx + 1
            end
        end
    end

    df = DataFrame(Dict("ka" => kas,
                        "ks" => kss,
                        "ki" => kis,
                        "ad" => ads,
                        "pk" => pks,
                        "pk2" => pk2s))
    if (! small) 
        CSV.write("../loop_simulations/20_response_peaks.csv", df)
    end
end

t_span = 3

function get_cascade(n_mediators, ks, ka, input_level)
    input = 1
    mediators = collect(2:1+n_mediators)
    control = n_mediators+2
    output = n_mediators+3

    sources = zeros(0)
    destinations = zeros(0)
    wghts = zeros(0)

    for m in mediators
        push!(sources, input)
        push!(destinations, m)
        push!(wghts, 1.0)

        push!(sources, m)
        push!(destinations, output)
        push!(wghts, -1.0)
    end

    if (length(mediators)==0)
        push!(sources, input)
        push!(destinations, output)
        push!(wghts, -1.0)
    end

    push!(sources, input)
    push!(destinations, control)
    push!(wghts, -1.0)

    sources = [Int(s) for s in sources]
    destinations = [Int(d) for d in destinations]

    digraph = SimpleWeightedDiGraph(sources, destinations, wghts)

    iid = 1
    spans = [[0, t_span]]

    u0 = zeros(0)
    append!(u0, [input_level, 1.0-input_level, 0.0])
    for m in mediators
        append!(u0, [0.0, 1.0, 0.0])
    end
    append!(u0, [1.0, 0.0, 0.0])
    append!(u0, [1.0, 0.0, 0.0])

    get_k_bkgd = get_k_bkgd_func(0.0)
    get_k_vals = get_kval_func(Dict('s' => ks, 'i' => 0), Dict('a' => ka, 'p' => 0))
    iff_parms = construct_parms(digraph, iid, spans, u0, Nothing, Nothing, get_k_bkgd, get_k_vals)
    sol = run_digraph(iff_parms, t_span, 0.001)
    return sol
end

function plot_solutions(solutions, labels, fname="none")
    @assert length(solutions) == length(labels)
  
    cscheme = get(ColorSchemes.viridis, reverse(0.1:0.8/(length(labels)):0.9))

    @. hfn(x,p) = p[1] + p[2] * ((x+0im)^p[3])/((x+0im)^p[3] + (p[4]+0im)^p[3])
    hcoeffs = zeros(length(labels))
    kd50s = zeros(length(labels))

    p = plot()
    for i in 1:length(labels)
        sol = solutions[i]
        lab = labels[i]
        xvals = sol[1]
        yvals = [y for y in sol[2][:, end]]

        fit = curve_fit(hfn, xvals, yvals, [1.0, -1.0, 1.0, 1.0])
        if fit.converged
            hcoeffs[i] = fit.param[3]
            kd50s[i] = fit.param[4]
        end

        if lab == "0 Repressor"
            p = plot!(p,
                      xvals, yvals,
                      label="Direct Repression", color=:grey,
                      linewidth=2, linestyle=:dash)
        else
            if lab != "1 Repressor"
                lab = lab * "s"
            end
            p = plot!(p,
                      xvals, yvals,
                      label=lab, color=cscheme[i],
                      linewidth=2, linestyle=:solid)
        end
    end

    p = plot!(p,
              #  color_palette=ColorSchemes.viridis,
              xlabel="Days of Recruitment",
              ylabel="Fraction Active",
              ylims=(-0.1, 1.1),
              yticks=0:0.5:1.0,
              xlims=(0, t_span),
              xticks=0:t_span/2.0:t_span,
              legend=:topright,
              #  legendtitle="# Repressors",
              size=(450, 300))

    if fname != "none"
        savefig(p, fname)
    end
    
    return [hcoeffs, kd50s]
end

function get_coeffs(n_mediators, ks, ka, input_level)
    sol = get_cascade(n_mediators, ks, ka, input_level)

    @. hfn(x,p) = 1 - ((x)^p[1])/((x)^p[1] + (p[2])^p[1])

    xvals = sol[1]
    yvals = [y for y in sol[2][:, end]]
    fit = curve_fit(
                    hfn, xvals, yvals, [1.0, 1.0], 
                    lower=[0.0, 0.0], autodiff=:forwarddiff)

    if fit.converged
        return fit.param
    else
        return [NaN, NaN]
    end
end

n_reps = [0, 2, 5, 20]
parms = [[n, 2, 1.0, 1.0] for n in n_reps]
cascades = [get_cascade(p...) for p in parms]
plot_solutions(cascades, [string(n)*" Repressor" for n in n_reps], "fig_4j.pdf");

function get_hill_dbase(small = false) 
    if ! small
        nmrange = [1, 2, 5, 10, 20]
        ksrange = 0:0.25:20
        karange = 0:0.1:2
        karange = 0:0.25:20
    else
        nmrange = 19:20
        ksrange = 19:20
        karange = 19:20
    end

    dflen = length(nmrange) * length(ksrange) * length(karange)
    println("Should take about " * string(dflen * 0.2) * " seconds")
    nms = zeros(dflen)
    kss = zeros(dflen)
    kas = zeros(dflen)
    ins = zeros(dflen)
    hcs = zeros(dflen)
    kds = zeros(dflen)
    idx = 1

    t1 = time()
    for n_m in reverse(nmrange)
        for ks in reverse(ksrange)
            for ka in reverse(karange)
                @time h, k = get_coeffs(n_m, ks, ka, 1.0)
                nms[idx] = n_m
                kss[idx] = ks
                kas[idx] = ka
                ins[idx] = 1.0
                hcs[idx] = h
                kds[idx] = k
                idx = idx + 1
            end
        end
        GC.gc()
    end
    t2 = time()

    #  println(t2 - t1)

    @assert idx - 1 == dflen

    cascade_df = DataFrame(Dict("nmediators" => nms,
                                "ks" => kss,
                                "ka" => kas,
                                "input" => ins,
                                "hillcoeff" => hcs,
                                "kd50" => kds))
    CSV.write("../loop_simulations/hill_coeffs.csv", cascade_df)
end

#  hill_dbase(true)
#  @time hill_dbase(false)


println("To generate full .csv files run:")
println("    get_response_dbase(false)")
println("    get_hill_dbase(false)")

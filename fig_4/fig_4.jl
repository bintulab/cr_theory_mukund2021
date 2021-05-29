using Colors, ColorSchemes
using DataFrames
using StatsPlots
using Statistics
using GLM
using Measures
using LsqFit

pyplot(grid=false)
fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)
#  Plots.resetfontsizes()

noise_level = 0.003544

###############################################################################
### FIG 5B
###############################################################################

# note n100 here really refers to 20-bin output discretization
#  hdac4_n100_ccs = Dict(0.0 => [0.8402011466907346, 2.0],
#                        10.0 => [0.005863478544182091, 0.030183266884945978],
#                        3.5 => [0.017704366746612615, 1.0086261681139732],
#                        12.5 => [0.0022227057791606012, 0.0],
#                        17.0 => [0.0015924353307981089, 0.327234862587496],
#                        5.0 => [0.007610191710702948, 1.134158217473606])
#
#  krab_n100_ccs = Dict(0.0 => [1.16067012194699, 1.952447944687942],
#                       10.0 => [0.0, 1.739110607715436],
#                       3.5 => [0.4876737948231341, 1.7115846964875094],
#                       12.5 => [0.1145708105550367, 1.7152622301665355],
#                       17.0 => [0.08268430173801551, 1.5],
#                       5.0 => [0.30676111081711127, 2.0])
#
hdac4_n100_ccs = Dict(0.0 => [[0.8395657535269593, 0.9398621506747077, 0.8698482155015244],
                              [2.0, 2.0, 2.0]],
                      10.0 => [[0.005875444366093918, 0.006125867332119918, 0.0058438457002036835],
                               [0.11183885603877655, 0.032062847509251396, 0.04804437896517706]],
                      3.5 => [[0.01861918152718385, 0.01366680253467901, 0.01724116224526335],
                              [1.0034538683709804, 1.1012132399802703, 1.05730693992646]],
                      12.5 => [[0.0025543496433731933, 0.0022316627076182944, 0.0018161094420840282],
                               [0.008057310112370875, 6.902924490377821e-5, 0.00028864160113290394]],
                      17.0 => [[0.001811719983757539, 0.001605769143032653, 0.0017818202331329332],
                               [0.33820056264903425, 0.3391229725623339, 0.36219747515886513]],
                      5.0 => [[0.008149041067623639, 0.007136221333670726, 0.008208155790609689],
                              [1.1243352519815528, 1.14223217198699, 1.0643099214235372]])
#
krab_n100_ccs = Dict(0.0 => [[1.1641148122813332, 1.1535267655391181, 1.1852105383519198],
                             [1.94655066828946, 1.8982960259224548, 1.8389020123437079]],
                     10.0 => [[0.16380093476272045, 0.16138959503175376, 0.16296417601343205],
                              [1.681785848519005, 1.6176462258163011, 1.7147928409538287]],
                     3.5 => [[0.5026229177748988, 0.4875715218096181, 0.48684338150625844],
                             [1.7241839285284164, 1.7321876361792694, 1.6386850451210508]],
                     12.5 => [[0.11302353845993598, 0.11666404031033872, 0.11452831927337837],
                              [1.7465270780278042, 1.5452832442279902, 1.6696732025861647]],
                     17.0 => [[0.08304292492697979, 0.08328258582350384, 0.08272769399505452],
                              [1.5, 1.5025220976978841, 1.5]],
                     5.0 => [[0.3083325885905498, 0.3057476879921178, 0.30717772830552587],
                             [2.0, 2.0, 2.0]])

krab_ccs = Dict(100 => krab_n100_ccs)
hdac4_ccs = Dict(100 => hdac4_n100_ccs)

krab_t_keys = collect(keys(krab_n100_ccs))
hdac_t_keys = collect(keys(hdac4_n100_ccs))

cr_list = ["HDAC4", "KRAB"]
cr_dict = Dict("HDAC4"  => hdac4_ccs,
               "KRAB"   => krab_ccs)
times_dict = Dict("HDAC4"   => hdac_t_keys,
                  "KRAB"   => krab_t_keys)
sample_types = ["Single-Cell", "Population"]
sample_idx_dict = Dict("Single-Cell" => 1,
                       "Population"  => 2)

function get_cc(cr, num_bins, sample_type, time)
    cc_t_dict = cr_dict[cr][num_bins][time][sample_idx_dict[sample_type]]
end

function get_plot_ccs(num_bins, sample_type, plot_size=(600, 600))
    krab_cc_dict   = krab_ccs[num_bins]
    hdac4_cc_dict  = hdac4_ccs[num_bins]

    if (sample_type=="Single-Cell")
        idx = 1
        upper_lim = 2
    else
        idx = 2
        upper_lim = 2
    end

    color_map = Dict("KRAB" => "darkorange",
                     "HDAC4" => "purple")

    x_hdac4 = sort(collect(keys(hdac4_n100_ccs)))
    y_hdac4 = [hdac4_cc_dict[t][idx] for t in x_hdac4]

    x_krab = sort(collect(keys(krab_n100_ccs)))
    y_krab = [krab_cc_dict[t][idx] for t in x_krab]

    # we fit an exponential model to these things
    model(t, p) = p[1] * exp.(-p[2] * t)
    kfit = curve_fit(model, x_krab, [mean(y) for y in y_krab], [0.5, 0.5])
    hfit = curve_fit(model, x_hdac4, [mean(y) for y in y_hdac4], [0.5, 0.5])

    kfunc(t) = kfit.param[1] * exp(-kfit.param[2] * t)
    hfunc(t) = hfit.param[1] * exp(-hfit.param[2] * t)

    ps = plot()
    ps = plot!(kfunc, x_krab[1], x_krab[end], 
               label = "KRAB",
               color = color_map["KRAB"],
               linewidth = 2)
    ps = plot!(hfunc, x_hdac4[1], x_hdac4[end],
               label = "HDAC4",
               color = color_map["HDAC4"],
               linewidth = 2)
    ps = plot!(x_krab, [mean(y) for y in y_krab], seriestype = :scatter,
               yerr = [std(y) for y in y_krab],
               label=false,
               color=color_map["KRAB"],
               markersize = 8,
               linewidth = 0.5)
    ps = plot!(x_hdac4, [mean(y) for y in y_hdac4], seriestype = :scatter,
               yerr = [std(y) for y in y_krab],
               label=false,
               color=color_map["HDAC4"],
               markersize = 8,
               linewidth = 0.5)

    #  for i in 1:length(y_krab[1])
    #      ps = plot!(x_krab, [y[i] for y in y_krab], seriestype=:scatter,
    #                 label="KRAB", markersize=6,
    #                 color=color_map["KRAB"])
    #      ps = plot!(ps,
    #                 x_hdac4, [y[i] for y in y_hdac4], seriestype=:scatter,
    #                 label="HDAC4", markersize=6,
    #                 color=color_map["HDAC4"],
    #                 xlabel="Days",
    #                 ylabel="Channel Capacity (bits)", ylims=(-0.1, upper_lim+0.1),
    #                 yticks=0:1.0:upper_lim,
    #                 legendfont = font(10),
    #                 size=plot_size)
    #  end
    ps = plot!(ps, 
               xlabel="Days Released",
               ylabel="Channel Capacity (bits)", ylims=(-0.1, upper_lim+0.1),
               yticks=0:1.0:upper_lim,
               legendfont = font(10),
               size=plot_size)
    return ps
end

pp = get_plot_ccs(100, "Population")
pp = plot!(pp, title="Population", top_margin=10mm)
ps = get_plot_ccs(100, "Single-Cell")
ps = plot!(ps, title="Single-Cell")

lay = @layout [a; b]
p = plot(ps, pp, layout=lay)
savefig(p, "fig_4c-d.pdf")

###############################################################################
### FIG 5C
###############################################################################

function get_pct_off(t_recruit, t_release, k)
    ks = k['s']
    ka = k['a']
    ki = k['i']
    y1 = (ks+ka+ki - sqrt((ks+ka+ki)^2 - 4*ks*ki))/2
    y2 = (ks+ka+ki + sqrt((ks+ka+ki)^2 - 4*ks*ki))/2

    function c_r(t_recruit)
        return ks/(y2-y1) * (exp(-y1*t_recruit) - exp(-y2*t_recruit))
    end

    function c_i(t_recruit)
        if (ki==0)
            return 0
        else
            term1 = (ks*ki)/(y1*(y2-y1)) * exp(-y1*t_recruit)
            term2 = (ks*ki)/(y2*(y2-y1)) * exp(-y2*t_recruit)
            return 1 - term1 + term2
        end
    end

    return c_r(t_recruit) * exp(-ka * t_release) + c_i(t_recruit)
end

# see https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1023694
function get_single_cell_cc(t_recruit, t_release, k)
    p = 1 - get_pct_off(t_recruit, t_release, k) # error probability
    return log2(1 + (1-p)^(p/(1-p)))
end

function get_population_cc(t_recruit, t_release, k, variance=0.05)
    power = get_pct_off(t_recruit, t_release, k)
    return maximum([0, 0.5 * log2(1 + power/variance)])
end

k_eed = Dict('a' => 0.05,
             's' => 0.9,
             'i' => 0.09)

k_hdac = Dict('a' => 0.76,
              's' => 8.6,
              'i' => 0)

k_krab = Dict('a' => 0.31,
              's' => 11,
              'i' => 0.13)

t_recruit = 5
t_releases = collect(0:0.01:15)

k_s_ccs = [get_single_cell_cc(t_recruit, t_release, k_krab) for t_release in t_releases]
k_p_ccs = [get_population_cc(t_recruit, t_release, k_krab, noise_level) for t_release in t_releases]
h_s_ccs = [get_single_cell_cc(t_recruit, t_release, k_hdac) for t_release in t_releases]
h_p_ccs = [get_population_cc(t_recruit, t_release, k_hdac, noise_level) for t_release in t_releases]

color_map = Dict("KRAB" => "darkorange",
                 "HDAC4" => "purple")

p1 = plot(t_releases, k_p_ccs,
          label="KRAB",
          color=color_map["KRAB"],
          linewidth=2)
p1 = plot!(t_releases, h_p_ccs,
           label="HDAC4",
           color=color_map["HDAC4"],
           linewidth=2,
           xlabel="Days Released",
           ylabel="Channel Capacity (bits)", 
           ylims=(-0.1, 4.1), yticks=0:1.0:4.0,
           title="Population",
           legend=:right,
           legendfont = font(10),
           top_margin=10mm,
           size=(600, 600))


p2 = plot(t_releases, k_s_ccs,
          label="KRAB",
          color=color_map["KRAB"],
          linewidth=2)
p2 = plot!(t_releases, h_s_ccs,
           label="HDAC4",
           color=color_map["HDAC4"],
           linewidth=2,
           xlabel="Days Released", 
           ylabel="Channel Capacity (bits)", ylims=(-0.1,2.1),
           yticks=0:1.0:3.0,
           legendfont = font(10),
           title="Single-Cell")
lay = @layout[a; b]
p = plot(p2, p1, layout=lay)
savefig(p, "fig_4g-h.pdf")

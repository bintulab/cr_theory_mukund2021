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
### FIG S3A
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

###############################################################################
### FIG 5D
###############################################################################


t_recruits = collect(0:0.001:1)
t_release = 0
k_p_ccs = [get_population_cc(t_recruit, t_release, k_krab, noise_level) for t_recruit in t_recruits]
h_p_ccs = [get_population_cc(t_recruit, t_release, k_hdac, noise_level) for t_recruit in t_recruits]

p1 = plot(t_recruits, k_p_ccs,
         label="KRAB",
         color=color_map["KRAB"],
         linewidth=2)
p1 = plot!(t_recruits, h_p_ccs,
          label="HDAC4",
          color=color_map["HDAC4"],
          linewidth=2,
          xlabel="Days Recruited", 
          ylabel="Channel Capacity (bits)", yticks=[0, 1, 2, 3, 4],
          xticks = [0, 0.5, 1],
          legend=:right,
          legendfont = font(10),
          size=(600, 300))

savefig(p1, "fig_s3a.pdf")


###############################################################################
### FIG 5E
###############################################################################


t_recruits = 5
t_release = 1
sigmas = collect(0.00001:0.00001:0.25)
k_p_ccs = [get_population_cc(t_recruit, t_release, k_krab, sigma) for sigma in sigmas]
h_p_ccs = [get_population_cc(t_recruit, t_release, k_hdac, sigma) for sigma in sigmas]

p1 = plot(sigmas, k_p_ccs,
         label="KRAB",
         color=color_map["KRAB"],
         linewidth=2)
p1 = plot!(sigmas, h_p_ccs,
          label="HDAC4",
          color=color_map["HDAC4"],
          linewidth=2,
          xlabel="Output Noise Variance", 
          # xscale=:log10,
          # legend=:outerright,
          ylabel="Channel Capacity (bits)",
          ylims=(0,8.5),
          legendfont = font(10),
          size=(600, 300))

println("KRAB: $(get_population_cc(t_recruit, t_release, k_krab, noise_level))")
println("HDAC4: $(get_population_cc(t_recruit, t_release, k_hdac, noise_level))")

p1 = vline!(p1, [noise_level], color=:black, linestyle=:dash, label=false)
savefig(p1, "fig_s3b.pdf")


###############################################################################
### FIG 5F
###############################################################################

function get_population_cc_input_noise(t_recruit, t_release, k, input_variance = 0.05, output_variance=0.05)
    target_poff = get_pct_off(t_recruit, t_release, k)

    # input b/w zero and 1 (power=1), and noise is input variance
    input_variance_scaled = target_poff * input_variance

    overall_cc = 0.5*log2(1 + target_poff/(input_variance_scaled + output_variance))

    return overall_cc
end

k_hdac = Dict('a' => 0.76,
              's' => 8.6,
              'i' => 0)

k_krab = Dict('a' => 0.31,
              's' => 11,
              'i' => 0.13)

t_recruit = 5
t_release = 5

sigmas = collect(0.00001:0.00001:0.25)

input_noise_sigmas = [0, 0.001, 0.01, 0.1, 1]
noise_colors = reverse([colorant"purple4", colorant"purple1", colorant"mediumorchid2", colorant"violet", colorant"plum1", colorant"thistle1"])
noise_colors = reverse(noise_colors)
noise_colors = get(ColorSchemes.viridis, reverse(range(0.0, stop=0.9, length=length(noise_colors))))

krab_cc_list = [[get_population_cc_input_noise(1, t_release, k_krab, input_noise_sigmas[i], sigma) for sigma in sigmas] for i in 1:length(input_noise_sigmas)]
hdac_cc_list = [[get_population_cc_input_noise(1, t_release, k_hdac, input_noise_sigmas[i], sigma) for sigma in sigmas] for i in 1:length(input_noise_sigmas)]

# t1_p_ccs = [get_population_cc_input_noise(1, t_release, k_krab, sigma) for sigma in sigmas]
# t2_p_ccs = [get_population_cc_input_noise(2, t_release, k_krab, sigma) for sigma in sigmas]
# t3_p_ccs = [get_population_cc_input_noise(3, t_release, k_krab, sigma) for sigma in sigmas]
# t4_p_ccs = [get_population_cc_input_noise(4, t_release, k_krab, sigma) for sigma in sigmas]
# t5_p_ccs = [get_population_cc_input_noise(5, t_release, k_krab, sigma) for sigma in sigmas]
#
# cc_list = [t1_p_ccs, t2_p_ccs]
p1 = plot(xlabel="Output Noise Variance", 
          ylabel="Channel Capacity (bits)",
          title="KRAB",
          legendtitle="Input Noise Variance",
          legend=false,
          legendfont = font(10),
          xscale=:log10,
          ylims=(-0.25,8.5),
          size=(600, 300))
for i in 1:length(krab_cc_list)
    p1 = plot!(sigmas, krab_cc_list[i], label="$(input_noise_sigmas[i])",
               color=noise_colors[i],
               linewidth=2)
end

p2 = plot(xlabel="Output Noise Variance", 
          ylabel="Channel Capacity (bits)",
          title="HDAC4",
          legendtitle="Input Noise Variance",
          legendfont = font(10),
          xscale=:log10,
          legend=:outerright,
          ylims=(-0.25,8.5),
          size=(1200, 300))
for i in 1:length(hdac_cc_list)
    p2 = plot!(sigmas, hdac_cc_list[i], label="$(input_noise_sigmas[i])",
               color=noise_colors[i],
               linewidth=2)
end

p1 = vline!(p1, [noise_level], color=:black, linestyle=:dash, label=false)
p2 = vline!(p2, [noise_level], color=:black, linestyle=:dash, label=false)

lay = @layout[a  b]
p = plot(p1, p2, layout=lay)
savefig(p, "fig_s3c.pdf")

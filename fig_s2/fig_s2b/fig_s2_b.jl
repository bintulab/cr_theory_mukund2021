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

fnt = Plots.font("Arial")
default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)

epsilon = 1e-8

################################################################################
####  FIG S2B -- write data to CSV
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

ka = 0.75
ks = 2
ki = 0.5

k_recruit = Dict('a' => ka,
                 's' => ks,
                 'i' => ki,
                 'p' => 0)
k_release = Dict('a' => ka,
                 's' => 0,
                 'i' => 0,
                 'p' => 0)

function get_a(rec, rel)
    ktp = [[k_recruit, (0.0, rec)],
           [k_release, (rec, rec+rel)]]
    sol = run_recruitment(ktp)
    a = sol[2][end, 1]
    return a
end

function get_i(rec, rel)
    ktp = [[k_recruit, (0.0, rec)],
           [k_release, (rec, rec+rel)]]
    sol = run_recruitment(ktp)
    a = sol[2][end, 3]
    return a
end

recs = zeros(0)
rels = zeros(0)
vals = zeros(0)
labs = String[]

for rec in 0.01:0.1:5.01
    for rel in 0.01:0.1:5.01
        #  println(string(rec)*" "*string(rel))
        push!(recs, rec)
        push!(rels, rel)
        push!(labs, "Active")
        push!(vals, get_a(rec, rel))

        push!(recs, rec)
        push!(rels, rel)
        push!(labs, "Irreversibly Silent")
        push!(vals, get_i(rec, rel))
    end
end

data = DataFrame(Dict("Recruitment Duration" => recs,
                      "Release Duration" => rels,
                      "Fraction" => vals,
                      "Population" => labs))

CSV.write("recruitment_graded_response.csv", data)

using Distributed
using DifferentialEquations
using ParameterizedFunctions
using Plots
using Random
using DataFrames
using CSV

pyplot()

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

function run_recruitment(k_time_pairs)
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
    data = [as rs is]

    return [t, data]
end

function get_sim_data(ka, ks, ki)
    kas = Float64[]
    kss = Float64[]
    kis = Float64[]
    kps = Float64[]
    ts = Float64[]
    as = Float64[]
    rs = Float64[]
    is = Float64[]
    k_recruit = Dict(
        'a' => ka,
        's' => ks,
        'i' => ki,
        'p' => 0)
    k_release = Dict(
        'a' => ka,
        's' => 0,
        'i' => 0,
        'p' => 0)
    recruit_pair = [k_recruit, (0.0, 5.0)]
    release_pair = [k_release, (5.0, 35.0)]
    sol = run_recruitment([recruit_pair, release_pair])
    ts = [ts ; sol[1]]
    as = [as ; sol[2][:, 1]]
    rs = [rs ; sol[2][:, 2]]
    is = [is ; sol[2][:, 3]]

    append!(kas, transpose([ka for i in range(1, stop=length(sol[1]))]))
    append!(kss, transpose([ks for i in range(1, stop=length(sol[1]))]))
    append!(kis, transpose([ki for i in range(1, stop=length(sol[1]))]))
    append!(kps, transpose([0  for i in range(1, stop=length(sol[1]))]))
    return [kas kss kis kps ts as rs is]
end

function run_sims_orig()
    k_scale_ki = collect(0:0.025:0.5)
    k_scale_ka = collect(0:0.05:1)
    k_scale_ks = collect(0:0.5:10)
    for ka in k_scale_ka
        for ks in k_scale_ks
            for ki in k_scale_ki
                data = get_sim_data(ka, ks, ki)
                df_dict = Dict(
                    "ka" => data[:, 1],
                    "ks" => data[:, 2],
                    "ki" => data[:, 3],
                    "kp" => data[:, 4],
                    "t" => data[:, 5],
                    "a" => data[:, 6],
                    "r" => data[:, 7],
                    "i" => data[:, 8])
                df = DataFrame(df_dict)
                println("./sims_orig_3sm_5d/ka-" * string(ka) * "_ks-" * string(ks) * "_ki-" * string(ki) *".csv")
                CSV.write("./sims_orig_3sm_5d/ka-" * string(ka) * "_ks-" * string(ks) * "_ki-" * string(ki) *".csv", df)
            end
        end
    end
end

function run_sims_rand(num_samples = 8000)
    function get_df(kvals)
        ka, ks, ki = kvals .* [1, 10, 0.5]
        data = get_sim_data(ka, ks, ki)
        df_dict = Dict(
                       "ka" => data[:, 1],
                       "ks" => data[:, 2],
                       "ki" => data[:, 3],
                       "kp" => data[:, 4],
                       "t" => data[:, 5],
                       "a" => data[:, 6],
                       "r" => data[:, 7],
                       "i" => data[:, 8])
        df = DataFrame(df_dict)
        return df
    end
    
    kvals = [rand(3) for i in 1:num_samples]
    df = vcat(pmap(get_df, kvals)...)
    CSV.write("./randomized_cr_df.csv", df)
    println("Wrote dataframe for " * string(num_samples) * " samples.")
end

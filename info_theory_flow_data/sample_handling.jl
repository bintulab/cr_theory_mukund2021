using Statistics

num_populations = 100

#= Given FITC array, returns percent off (i.e. < 10^4)
=#
function get_pct_off(fitcs)
    count = 0

    if length(fitcs)==0
        println("no samples!")
    end

    for i in fitcs
        if i < 10000
            count = count + 1
        end
    end

    pct = count/length(fitcs) * 100.0

    #  if pct > 1
    #      println(pct)
    #      println(length(fitcs))
    #      println(count)
    #      println()
    #  end
    #
    return pct
end

#= Given dict of input => FITCs, returns a dict of input => pct off =#
function get_pct_off_subsample(fitc_vals_dict)
    return Dict(i => get_pct_off(fitc_vals_dict[i])
                for i in collect(keys(fitc_vals_dict)))
end

#= Given array of samples and array of keys, splits sample by keys

Example: arr=[1,2,3,4], ks=[1,2,1,2] ==> [[1,3], [2,4]]
=#
function stratify_sample(arr, ks)
    d = Dict(i => Float64[] for i in sort(unique(ks)))
    for i in 1:length(ks)
        append!(d[ks[i]], arr[i])
    end

    return [d[i] for i in sort(collect(keys(d)))]
end

#= Splits a population (dict key => all_members) into subpopulations

pop_dict is a dictionary from keys => all members corresponding to that key
num_populations is number of subpopulations to split the population into

Returns dict of subpop_idx => key => subpop_members
=#
function subsample_population(pop_dict, num_populations)
    function stratify_entry(pop_dict_entry, num_pops)
        subpop_keys = rand(1:num_pops, length(pop_dict_entry)) #ensure partitioning
        if length(unique(subpop_keys)) < num_populations # make sure each subpop gets members
            while (length(unique(subpop_keys)) < num_populations)
                subpop_keys = rand(1:num_pops, length(pop_dict_entry))
            end
        end
        return stratify_sample(pop_dict_entry, subpop_keys)
    end

    # give each member for each key a subpopulation index
    pop_split_dict = Dict(i => stratify_entry(pop_dict[i], num_populations)
                            for i in collect(keys(pop_dict)))

    # build subpopulation dictionary
    pop_idx = collect(range(1, stop=num_populations))
    return map(idx -> Dict(i => pop_split_dict[i][idx]
                            for i in collect(keys(pop_split_dict))),
                pop_idx)
end


function subsample_population_bioreps(pop_dict)
    # given
    #   input ==> biorep ==> [samples]
    # want to return
    #   biorep ==> input ==> [samples]

    samples = pop_dict
    inputs = collect(keys(samples))
    bioreps = collect(keys(samples[inputs[1]]))

    return Dict(biorep =>
                Dict(input => samples[input][biorep]
                     for input in inputs)
                for biorep in bioreps)
end

function partition_bioreps(biorep_dict)
    # given 
    #   input ==> biorep ==> samples
    # want to return
    #   idx ==> input ==> samples
    # where each biorep has been split into 10 idx's 
    samples_per_biorep = 30
    inputs = collect(keys(biorep_dict))
    bioreps = collect(keys(biorep_dict[inputs[1]]))

    pdict = Dict(idx => 
                 Dict(input => Float64[] for input in inputs)
                 for idx in 1:samples_per_biorep*length(bioreps))

    for biorep in bioreps
        for input in inputs
            b_samples = biorep_dict[input][biorep]

            strat_keys = rand(1:samples_per_biorep, length(b_samples))
            if length(unique(strat_keys)) < samples_per_biorep # make sure each subpop gets members
                while (length(unique(strat_keys)) < samples_per_biorep)
                    strat_keys = rand(1:samples_per_biorep, length(b_samples))
                end
            end

            stratified_bsamp = stratify_sample(b_samples, strat_keys)

            idx_base = samples_per_biorep * (biorep-1)

            for i in 1:samples_per_biorep
                pdict[idx_base + i][input] = stratified_bsamp[i]
            end
        end
    end

    return pdict
end

function bootstrap_bioreps(biorep_dict)
    # given 
    #   input ==> biorep ==> samples
    # want to return
    #   idx ==> input ==> samples
    # where each biorep has been bootstrapped into 100 subsamples
    
    num_bootstrapped_samples = 100
    inputs = collect(keys(biorep_dict))
    bioreps = collect(keys(biorep_dict[inputs[1]]))

    pdict = Dict(idx => 
                 Dict(input => Float64[] for input in inputs)
                 for idx in 1:num_bootstrapped_samples*length(bioreps))

    for biorep in bioreps
        for input in inputs
            idx_base = num_bootstrapped_samples*(biorep-1)
            b_samples = biorep_dict[input][biorep]
            bootstraps = [sample(b_samples, length(b_samples), replace=true) for i in 1:num_bootstrapped_samples]
            for i in 1:num_bootstrapped_samples
                pdict[idx_base + i][input] = bootstraps[i]
            end
        end
    end
    return pdict
end

#= Given dict of input => biorep => FITCs, returns population level channel capacity

Summary statistic used is % off
=#
function get_pop_cc(date_dict)
    # input => biorep => FITC[]
    #samples = get_samples(date_dict)
    samples = date_dict

    # subpops = subsample_population(samples, num_populations[>num_populations<]) #partition into 100 groups

    subpops = subsample_population_bioreps(samples)

    pct_offs = Dict(i => get_pct_off_subsample(subpops[i])
                    for i in collect(keys(subpops)))

    # dict of subpop index -> input -> pct off
    pct_off_dict = Dict(
        time => [pct_offs[i][time] for i in collect(keys(pct_offs))]
        for time in collect(keys(pct_offs[collect(keys(pct_offs))[1]])))
   
    #  println(pct_off_dict)

    # now we can split the pct off dict to get channel capacities 
    # subsample_dict = partition_bioreps(samples)
    subsample_dict = bootstrap_bioreps(samples)
    pct_offs = Dict(i => get_pct_off_subsample(subsample_dict[i])
                    for i in collect(keys(subsample_dict)))
    pct_off_dict = Dict(
        time => [pct_offs[i][time] for i in collect(keys(pct_offs))]
        for time in collect(keys(pct_offs[collect(keys(pct_offs))[1]])))

    var_dict = Dict(t => var([p/100.0 for p in pct_off_dict[t]]) for t in keys(pct_off_dict))
    var_dict = Dict(t => mean(var_dict[t]) for t in keys(var_dict))
    var_mean = mean([var_dict[t] for t in keys(var_dict)])
    println("Mean variance: $(var_mean)")

    # println(pct_off_dict)
    pop_ccs = compute_channel_capacity(pct_off_dict, true)

    return pop_ccs
end

#= Given dict of input => biorep => FITCs, returns single-cell level channel capacity

Takes the log of single-cells so binning is uniform in log space
=#
function get_sc_cc(date_dict)
    #samples = get_samples(date_dict)

    # pool samples together: go from inpput => biorep = [sample] to input => sample
    function pool_samples(biorep_dict)
        arr = []
        for k in collect(keys(biorep_dict))
            append!(arr, biorep_dict[k])
        end
        return arr
    end

    samples = [pool_samples(date_dict[biorep]) for biorep in collect(keys(date_dict))]

    samples_logged = Dict(
        i => map(log10, filter(s->s>0, samples[i]))
            for i in collect(keys(samples)))
    return compute_channel_capacity(samples_logged, false)
end

#=
# Turn dict of input => biorep => FITCs to subsample => input => FITCs
=#
function get_split_sc_cc(date_dict)
    # input => biorep => FITCs to subsample => input => biorep => FITCs 
    split_date_dict = split_pop_date_dict(date_dict)

    # pool samples together: turn a dict or index => [arr] to [array]
    function pool_samples(input_dict)
        arr = []
        for k in collect(keys(input_dict))
            append!(arr, input_dict[k])
        end
        return arr
    end

    split_indices = collect(keys(split_date_dict))
    inputs = collect(keys(split_date_dict[split_indices[1]]))
    split_date_dict = Dict(s => Dict(i => pool_samples(split_date_dict[s][i])
                                     for i in inputs)
                           for s in split_indices)
    log_date_dict = Dict(s => Dict(i => map(log10, filter(s->s>0, split_date_dict[s][i]))
                                   for i in inputs)
                         for s in split_indices)

    ccs = []

    for s in split_indices
        append!(ccs, compute_channel_capacity(log_date_dict[s], false))
    end
   
    return ccs
end

#= Given dict of input => biorep => FITCs, splits up the FITCs such that we have
# 3 separate dictionaries of input => biorep => FITCs
# want to return subsample => input => biorep => FITCs
=#
function split_pop_date_dict(date_dict)
    # get input => subsample => biorep => FITCs 
    split_date_dict = Dict(i => subsample_population(date_dict[i], 3) for i in collect(keys(date_dict)))
    inputs = collect(keys(date_dict))
    bioreps = collect(keys(date_dict[inputs[1]]))
    subsamples = [1, 2, 3]

    split_date_dict = Dict(s => Dict(i => Dict(b => split_date_dict[i][s][b] 
                                               for b in bioreps)
                                     for i in inputs)
                           for s in subsamples)
    return split_date_dict
end

function get_split_pop_cc(date_dict)
    split_dict = split_pop_date_dict(date_dict)
    split_indices = collect(keys(split_dict))
    ccs = []
    for i in split_indices
        append!(ccs, get_pop_cc(split_dict[i]))
    end
    return ccs
end

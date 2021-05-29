using Colors
using CSV
using DataFrames
using Distributed
using Distributions
using GLM
using InformationMeasures
using MultivariateStats
using Plots
using Statistics

# num_bins is in infotheory_handling.jl
# num_populations is in sample_handling.jl

debug = false

include("fcsfile_handling.jl")
include("infotheory_handling.jl")
include("sample_handling.jl")

tubelist_location = "./tubelist_wiggled_br1-2.csv" # use wiggled tube list
println("Reading in tubelist...")
tubelist = CSV.read(tubelist_location, DataFrame)

crs = unique(tubelist[!, Symbol("Regulator")])
bioreps = unique(tubelist[!, Symbol("Biorep")])
inputs = [0, 2, 3, 5]

# timepoint_dict = Dict(
#     "KRAB" => [0, 1, 2, 5, 7, 10, 13, 24, 31],
#     "EED" => [0, 1, 2, 3, 7, 16, 21, 24, 27, 31],
#     "DNMT3B" => [0, 1, 2, 4, 5, 7, 10, 28, 31],
#     "HDAC4" => [0, 1, 2, 5, 7, 10, 13, 24, 31])

timepoint_dict = Dict(
    "KRAB" => [0, 3.5, 5, 10, 12.5, 17],
    "EED" => [0, 3, 9.5, 12, 15.5],
    "DNMT3B" => [0, 3.5, 5, 10, 13, 17.5],
    "HDAC4" => [0, 3.5, 5, 10, 12.5, 17])


cr_input_dict = Dict(cr => 
                     Dict(timepoint => 
                          Dict(input => 
                               Dict(biorep => 
                                    String[] for biorep in bioreps) 
                               for input in inputs) 
                          for timepoint in timepoint_dict[cr]) 
                     for cr in crs)

#= Returns a dictionary of input => output fcs files from given tubelist csv

INPUTS: tubelist_location is a string that gives the path to a csv file containing
a list of tubes/conditions

OUTPUTS: cr_input_dict is a dict from cr => timepoint => input => [fcsfiles]
len(fcsfiles) in that is just the # of bioreps
=#
function get_cr_input_dict()
    for cr in crs
        for timepoint in timepoint_dict[cr]
            for input in inputs
                for biorep in bioreps
                    tubedf = tubelist[
                                    (tubelist[!, Symbol("Regulator")].==cr).&
                                    (tubelist[!, Symbol("Biorep")].==biorep).&
                                    (tubelist[!, Symbol("Days after dox")].==timepoint).&
                                    (tubelist[!, Symbol("Days with dox")].==input),
                                : ]
                    if length(tubedf[!, Symbol("Tube")]) != 1
                        println("WHOOPS: ")
                        print(cr)
                        print(", biorep ")
                        print(biorep)
                        print(", +dox ")
                        print(input)
                        print(", -dox ")
                        print(timepoint)
                        println(tubedf[!, Symbol("Tube")])
                    end
                    fcsfile = "./fcs/"*tubedf[!, Symbol("Tube")][1]
                    if fcsfile[length(fcsfile)-3:end]!=".fcs"
                        fcsfile = fcsfile*".fcs"
                    end
                    cr_input_dict[cr][timepoint][input][biorep] = [fcsfile]
                end
            end
        end
    end

    return cr_input_dict
end

println("Getting CR input dict...")
cr_input_dict = get_cr_input_dict()
println("Done.")

#= Given dict of input => fcs paths, returns dict of input => fitc vals =#
function get_input_fitc_dict(input_fcs_dict)
    # assumes input_fcs_dict is a dict from input => fcs paths

    input_biorep_fitc_dict = Dict(input=>
                                  Dict(biorep=>
                                       Float64[] for biorep in keys(input_fcs_dict[input]))
                                  for input in keys(input_fcs_dict))

    # input_fitc_dict = Dict(k=>Float64[] for k in keys(input_fcs_dict))

    for input in keys(input_fcs_dict)
        for (biorep, fcsfile) in input_fcs_dict[input]
            @assert length(fcsfile)==1
            input_biorep_fitc_dict[input][biorep] = get_filtered_flow_data(fcsfile[1])
        end
    end

    # for (input, fcsfilelist) in input_fcs_dict
    #     for fcsfile in fcsfilelist
    #         fitcs = get_filtered_flow_data(fcsfile)
    #         append!(input_fitc_dict[input], fitcs)
    #     end
    # end

    return input_biorep_fitc_dict
end
#= Given input CR dict, returns the mapping to FITCs with separated bioreps =#
function get_cr_fitc_dict(cr_input_dict)
    cr_fitc_dict = Dict(cr => 
                        Dict(timepoint => 
                             Dict(input => 
                                  Dict(biorep => Float64[] for biorep in bioreps)
                                  for input in inputs) 
                             for timepoint in timepoint_dict[cr]) 
                        for cr in crs)

    for cr in crs
        for timepoint in timepoint_dict[cr]
            cr_fitc_dict[cr][timepoint] = get_input_fitc_dict(cr_input_dict[cr][timepoint])
        end
    end

    return cr_fitc_dict
end

# function get_channel_capacities(cr_fitc_dict)
#     cr_cc_dict = Dict(cr =>
#                       Dict(timepoint => Float64[] for timepoint in timepoint_dict[cr])
#                       for cr in crs)
#
#     for cr in crs
#         for timepoint in timepoint_dict[cr]
#             fitc_dict = cr_fitc_dict[cr][timepoint]
#             pop_cc = get_pop_cc(fitc_dict)
#             snc_cc = get_sc_c(fitc_dict)
#             cr_cc_dict[cr][timepoint] = [snc_cc, pop_cc]
#         end
#     end
#
#     return cr_cc_dict
# end

#= Gets dict from timepoint => input => biorep => FITC[] from input_dict, cr =#
function get_fitcs_for_cr(cr_input_dict, cr)
    fitc_dict = Dict(timepoint => 
                     Dict(input => 
                          Dict(biorep => Float64[] for biorep in bioreps) 
                          for input in inputs) 
                     for timepoint in timepoint_dict[cr])

    for timepoint in timepoint_dict[cr]
        fitc_dict[timepoint] = get_input_fitc_dict(cr_input_dict[cr][timepoint])
    end

    return fitc_dict
end


function get_ccs_for_timepoints(timepoint_input_fitc_dict)
    timepoints = keys(timepoint_input_fitc_dict)
    cc_dict = Dict(timepoint => [] for timepoint in timepoints)
    for timepoint in timepoints
        input_dict = timepoint_input_fitc_dict[timepoint]
        sc_cc = get_split_sc_cc(input_dict)
        print(timepoint)
        println(" => ")
        pop_cc = get_split_pop_cc(input_dict)
        cc_dict[timepoint] = [Float64.(sc_cc), Float64.(pop_cc)]
    end
    return cc_dict
end

function get_ccs_for_cr(cr)
    fitcs = get_fitcs_for_cr(cr_input_dict, cr)
    ccs = get_ccs_for_timepoints(fitcs)
    return ccs
end

println("HDAC4: ")
println(get_ccs_for_cr("HDAC4"))
println()
println()
println("KRAB: ")
println(get_ccs_for_cr("KRAB"))

using Distributions
using LinearAlgebra
using StatsBase

num_bins = 20

function get_hist_edges(input_dict, nbins, population)
    all_outputs = reduce(vcat, [input_dict[i] for i in keys(input_dict)])

    max_output = maximum(all_outputs)
    min_output = minimum(all_outputs)
    
    if (population) # percent off
        # println("Population CC bin edges!")
        bin_edges = collect(-1:102.0/nbins:101.0)
    else # FITC
        # println("Single cell CC bin edges!")
        bin_edges = collect(min_output-1:(max_output+2-min_output)/nbins:max_output+1)
    end

    h = fit(Histogram, all_outputs, bin_edges)
    return collect(h.edges[1])
end

function find_bin_index(value, bin_edges) # wheeee no bounds checking here
    if debug
        print("Bin edges: ")
        println(bin_edges)
        print("Value: ")
        println(value)
    end
    for i in 1:length(bin_edges)-1
        if value >= bin_edges[i] && value < bin_edges[i+1]
            if debug
                print("Found index at: ")
                println(i)
            end
            return i
        else
            if debug
                print("Did not find index at ")
                println(i)
            end
        end
    end
end

# each row of matrix is input, each column is output
function build_conditional_matrix(input_dict, population)
    inputs = collect(keys(input_dict))
    output_bin_edges = get_hist_edges(input_dict, num_bins, population) # lower edges, bin into 20 things
    freqs = zeros(length(inputs), length(output_bin_edges)-1)

    if debug
        println("Frequencies: ")
        println(freqs)
    end

    for i in 1:length(inputs)
        input = inputs[i]
        for j in 1:length(input_dict[input])
            output = input_dict[input][j]
            bin_idx = find_bin_index(output, output_bin_edges)

            if isnothing(bin_idx)
                println(output_bin_edges)
                println(output)
            end
            if debug
                print("Adding count to ")
                print(i)
                print(", ")
                println(bin_idx)
            end
            freqs[i, bin_idx] += 1
        end
    end

    if debug
        println("Frequencies: ")
        println(freqs)
    end

    # now we have to normalize each row
    for i in 1:size(freqs)[1]
        tot = sum(freqs[i,:])
        for j in 1:length(output_bin_edges)-1
            freqs[i, j] = freqs[i, j]/tot
        end
    end

    if debug
        println("Frequencies: ")
        println(freqs)
    end

    return freqs
end

function blahut_arimoto(x)
    numrows, numcols = size(x)

    # check for negative entries
    for i in 1:numrows
        for j in 1:numcols
            if (x[i,j]<0)
                println("ERROR: negative entry in input matrix")
                return
            end
        end
    end

    # check for zero columns
    for i in 1:numcols
        arr = x[:,i]
        if arr == zeros(1, numrows)
            println("ERROR: there is a zero column in the matrix")
            return;
        end
    end

    # check for zero rows and make sure they sum to 1
    for i in 1:numrows
        arr = x[i,:]
        if arr == zeros(1, numcols)
            println("ERROR: there is a zero row in the matrix")
            return;
        end
        x[i,:] = x[i,:]/sum(x[i,:]) # normalize rows
    end

    r = ones(numrows)/numrows
    q = zeros(numrows, numcols)
    error_tolerance = 1e-5/numrows
    r1 = zeros(numrows)

    for i in 1:numrows
        x[i,:] = x[i,:]/sum(x[i,:])
    end
    for iter in 1:10000
        for j in 1:numcols
            q[:,j] = r .* x[:,j]
            q[:,j] = q[:,j]/sum(q[:,j])
        end
        for i in 1:numrows
            r1[i] = prod(q[i,:] .^ x[i,:])
        end
        r1 = r1/sum(r1)
        r = r1

        #  if (norm(r1-r) < error_tolerance)
        #      break
        #  else
        #      r = r1
        #  end
    end

    C = 0
    for i in 1:numrows
        for j in 1:numcols
            if r[i] > 0 && q[i,j] > 0
                C = C + r[i]*x[i,j]*log(q[i,j]/r[i])
            end
        end
    end

    C = C / log(2)
    return (C, r)
end

#= Computes the channel capacity given a dict of input => outputs
=#
function compute_channel_capacity(input_dict, population=true)
    m = build_conditional_matrix(input_dict, population)
    C, r = blahut_arimoto(m)
    return C
end

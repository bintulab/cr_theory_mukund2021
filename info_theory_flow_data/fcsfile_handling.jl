using FileIO
using FCSFiles
using GeometricalPredicates

#= Gets flow data for a dictionary of tubes.

Given a dictionary of index => fcs file name, returns a dictionary of
index => filtered flow values.
=#
function get_facs_data(tubedict)
    # want to use B1-A for information -- FITC
    facs_data = Dict(i => get_filtered_flow_data(tubedict[i]) for i in collect(keys(tubedict)))

    if debug
        println("got facs data")
    end

    return facs_data
end

#= Filters flow data for live cells

Takes arguments for FSC-A, SSC-A, and B1A -- uses a polygon gate made in
Cytoflow to filter B1A values if the corresponding FSC-A and SSC-A values
indicate that the cell is a live singlet.
=#
function filter_flow_data(fscs, sscs, mchs, b1a)
    min_coords = [minimum(fscs), minimum(sscs)]
    max_coords = [maximum(fscs), minimum(sscs)]

    # GeometricalPredicates requires all coordinates to be in [1, 2)
    fscs_scaled = [(f-min_coords[1])/(max_coords[1]-min_coords[1]) + 1.0 for f in fscs]
    sscs_scaled = [(s-min_coords[2])/(max_coords[2]-min_coords[2]) + 1.0 for s in sscs]

    point1 = (50000,  10000)
    point2 = (200000, 10000)
    point3 = (200000, 100000)
    point4 = (50000,  100000)
    point5 = (50000,  10000)

    #  point1 = (113097.94767349683, 188528.74478919798)
    #  point2 = (158699.69036009547, 77818.35897898073)
    #  point3 = (86319.3037719132,   32617.497996363905)
    #  point4 = (63789.55549038625,  44560.87239599226)
    #  point5 = (63789.55549038625,  104693.81475858176)
    #  point6 = (90599.38825827425,  179128.20042318373)
    #  point7 = (113097.94767349683, 188528.74478919798) # close the polygon
    #  points = [point1, point2, point3, point4, point5, point6]

    # Scales a point to be in the domain [1, 2) for use with GeometricalPredicates
    function scale_point(point, min_coords, max_coord)
        return GeometricalPredicates.Point(
                    (point[1]-min_coords[1])/(max_coords[1]-min_coords[1]) + 1.0,
                    (point[2]-min_coords[2])/(max_coords[2]-min_coords[2]) + 1.0)
    end

    point1 = scale_point(point1, min_coords, max_coords)
    point2 = scale_point(point2, min_coords, max_coords)
    point3 = scale_point(point3, min_coords, max_coords)
    point4 = scale_point(point4, min_coords, max_coords)
    point5 = scale_point(point5, min_coords, max_coords)
    #  point6 = scale_point(point6, min_coords, max_coords)
    #  point7 = scale_point(point7, min_coords, max_coords)

    livegate = Polygon(point1, point2, point3, point4, point5)
    #  livegate = Polygon(point1, point2, point3, point4, point5, point6, point7)

    b1a_filtered = Float64[]

    for i in 1:length(b1a)
        point = GeometricalPredicates.Point(fscs_scaled[i], sscs_scaled[i])
        if inpolygon(livegate, point)>0
            if (mchs[i] > 1000) # gate on mCherry-positive cells
                append!(b1a_filtered, b1a[i])
            end
        end
    end

    pct = (1.0 * length(b1a_filtered))/(1.0 * length(fscs_scaled))

    return b1a_filtered
end

#= Returns filtered flow data given path to FCS file

Accepts string path to FCS file, and returns the filtered flow data -- in this
case, we're only really concerned with B1-A values gated on FSC-A/SSC-A
=#
function get_filtered_flow_data(fcsfile)
    filedata = load(fcsfile)
    b1a = filedata["B1-A"]
    fsc = filedata["FSC-A"]
    ssc = filedata["SSC-A"]
    mchs = filedata["Y2-A"]

    return filter_flow_data(fsc, ssc, mchs, b1a)
end

#= Given a dictionary of input => FCS filepath, returns an array of [inputs. outputs]

Note here that each output is a specific value/read from the FCS file
=#
function get_data_array(facsdict)
    data_arr = []
    for input in collect(keys(facsdict))
        yvals = facsdict[input]
        append!(data_arr, [[input y] for y in yvals])
    end

    return [[[Float32(x[1]) for x in data_arr]]
            [[Float32(x[2]) for x in data_arr]]]
end

#= Given a dictionary of input => FCS files, returns a dict of input => samples =#
function get_samples(tubedict)
    samples_dict = Dict(
        i => get_filtered_flow_data(tubedict[i])
                    for i in collect(keys(tubedict)))
end

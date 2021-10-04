using GeoData, ArchGDAL, Dates, Pkg

function load_data(basedir, carrycap)
    # The directory for this project
    basedir = dirname(Pkg.project().path)
    # Load growthrates from tif files saved in growthmaps.jl
    datapath = joinpath(basedir, "data", "growthrates")
    # Set the timespan of the data. This will be looped by DynamicGrids
    # for years before or after this period.
    growthtimespan = DateTime(2020, 1):Month(1):DateTime(2020, 12)
    # Load, cut out Australia, and reorient the growth rate data
    aus = X(Between(113.0, 154.0)), Y(Between(-44.0, -10.0))
    rate_h, rate_p = map((:host, :parasitoid)) do species
        slices = map(readdir(joinpath(datapath, string(species)); join=true)) do path
            GeoArray(path; mappedcrs=EPSG(4326))[aus..., Band(1)] |>
                a -> replace_missing(a, 0.0f0) |>
                a -> permutedims(a, (Y, X))
        end
        A = cat(slices...; dims=Ti(growthtimespan))
    end

    # Take an array slice to use as a template
    A = rate_h[Ti(1)]

    # Make a masking layer for ocean areas
    mask = boolmask(A)

    # Create an initial population array for the host
    init_h = zero(A)
    # Add a host population around Cairns in Queensland
    cairns = Y(Between(-17, -20)), X(Between(144.0, 145.0))
    init_h[cairns...] .= carrycap

    # Initialise scattered parasitoid populations
    init_p = (_ -> (r = rand(); r < 0.1 ? carrycap * 1f-4 : 0.0f0)).(A)

    return (; init_h, init_p, rate_h, rate_p, mask)
end

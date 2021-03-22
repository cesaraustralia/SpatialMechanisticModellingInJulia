using GeoData, ArchGDAL, Dates

# The directory for this project
basedir = realpath(joinpath(@__FILE__, ".."))
# Load growthrates from tif files saved in growthmaps.jl
datapath = joinpath(basedir, "output", "growthrates")
# Set the timespan of the data. This will be looped by DynamicGrids
# for years before or after this period.
growthtimespan = DateTime(2020, 1):Month(1):DateTime(2020, 12)
# Load, cut out Australia, and reorient the growth rate data
aus = Lon(Between(113.0, 154.0)), Lat(Between(-44.0, -10.0))
rate_h, rate_p = map((:host, :parasite)) do species
    slices = map(readdir(joinpath(datapath, string(species)); join=true)) do path
        GDALarray(path; mappedcrs=EPSG(4326))[aus..., Band(1)] |>
            a -> replace_missing(a, 0.0f0) |>
            a -> permutedims(a, (Lat, Lon))
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
cairns = Lat(Between(-17, -20)), Lon(Between(144.0, 145.0))
init_h[cairns...] .= carrycap

# Initialise scattered parasite populations
init_p = (_ -> (r = rand(Float32); r < 0.1 ? carrycap * 1f-4 : 0.0f0)).(A)

# Define initialisation and auxilary data for the simulations.
# :rand is only required on GPU
initdata = (; H=init_h, P=init_p, rand=rand(Float32, size(init_h)))
auxdata = (; rH=rate_h, rP=rate_p)

# Set the time span
tspan = DateTime(2020, 1):Week(1):DateTime(2023, 1)

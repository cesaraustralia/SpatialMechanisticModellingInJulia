using GeoData, ArchGDAL, Dates

# Load growthrates from tif files
path = "/home/raf/.julia/dev/MethodsPaper/output/growthrates/"
# Set the timespan of the data. This will be looped by DynamicGrids
# for years before or after this period.
growthtimespan = DateTime(2020, 1):Month(1):DateTime(2020, 12)
aus = Lon(Between(113.0, 154.0)), Lat(Between(-44.0, -10.0))
using GeoData, ArchGDAL, Dates

rH, rP = map((:host, :parasite)) do species
    slices = map(readdir(joinpath(path, string(species)); join=true)) do path
        GDALarray(path; mappedcrs=EPSG(4326))[aus..., Band(1)] |>
            a -> replace_missing(a, 0.0f0) |>
            a -> permutedims(a, (Lat, Lon))
    end
    A = cat(slices...; dims=Ti(growthtimespan))
end

slice = rH[Ti(1)]
tspan = DateTime(2020, 1):Day(7):DateTime(2023, 1)

# A masking layer for ocean areas
mask = boolmask(slice)

# Add a host population in cairns
hostinit = zero(slice)
cairns = Lat(Between(-17, -20)), Lon(Between(144.0, 145.0))
hostinit[cairns...] .= carrycap

# A scattered parasite populations everywhere
parainit = (_ -> (r = rand(Float32); r < 0.1 ? carrycap * 1f-4 : 0.0f0)).(slice)

# Define initialisation and auxilary data for the simulations.
# :rand is only required on GPU
initdata = (; H=hostinit, P=parainit, rand=rand(Float32, size(hostinit)))
auxdata = (; rH=rH, rP=rP)

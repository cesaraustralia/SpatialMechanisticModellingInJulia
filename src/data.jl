using GeoData, ArchGDAL, Dates

# Load growthrates from tif files
path= "/home/raf/.julia/dev/MethodsPaper/output/growthrates/"
growthtimespan = DateTime(2017,1):Month(1):DateTime(2017,12)
aus = Lon(Between(113.3402, 153.9523)), Lat(Between(-43.62234, -10.65125))
rH, rP = map((:host, :parasite)) do species
    slices = map(readdir(joinpath(path, string(species)); join=true)) do path
        GDALarray(path; mappedcrs=EPSG(4326))[aus..., Band(1)] |>
            a -> replace_missing(a, 0.0f0) |>
            a -> permutedims(a, (Lat, Lon))
    end
    A = cat(slices...; dims=Ti(growthtimespan))
end

tspan = DateTime(2020, 1):Day(7):DateTime(2035, 1)
mask = boolmask(rH[Ti(1)])

hostpop = zero(rH[Ti(1)])
cairns = Lat(Between(-17, -20)), Lon(Between(144.0, 145.0))
hostpop[cairns...] .= carrycap

parapop = (a -> (r = rand(Float32); r < 0.1 ? carrycap : 0.0f0)).(hostpop)

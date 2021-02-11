#### Simple Growth Model #####################

using GrowthMaps, ArchGDAL, GeoData, RasterDataSources, Dates, Unitful, Plots

# Need to set a path for downloaded raster files
# ENV["RASTERDATASOURCES_PATH"] = "/home/YOUR_USERNAME/Data/"

# Define julia object to hold parameters to for the simple growth model 
# a struct with parameters, that inherits from GrowthModel
struct IntrinsicGrowth{Ic,U,If,D} <: GrowthModel
    intercept::Ic
    up::U
    inflection::If
    down::D
end

# Add a Growthmaps.rate method to calculate growth for a specific temperature
# Weird: we are adding a method to a function in the GrowthMaps mackage
# for the struct we have defined above.
function GrowthMaps.rate(m::IntrinsicGrowth, temp) 
    growthrate(m.intercept, m.up, m.inflection, m.down, temp)
end

growthrate(s, a, b, c, x) = x < c ? s + a*x/c - a : s + b - b*x/c  

# Define out growth rate models:
host_growth = IntrinsicGrowth(0.3, 9.0, 50.0, 35u"°C" |> u"K")
para_growth = IntrinsicGrowth(0.3, 9.0, 40.0, 25u"°C" |> u"K")

# Plot growth rate curves:
temps = collect(-15.0u"°C":0.1u"K":40.0u"°C")
p = plot(x -> GrowthMaps.rate(host_growth, x), temps; label="host")
plot!(p, x -> GrowthMaps.rate(para_growth, x), temps; label="parasite")
display(p)

# Specify the WorldClim Climate layer `tavg` with units in degrees C.
# In practise we would use a finer resolution dataset to capture the 
# range of temperature stresses, but this is easier to download and use.
host_model = Model(Layer(:tavg, u"°C", host_growth))
para_model = Model(Layer(:tavg, u"°C", para_growth))

# Download and set up the WorldClim Climate data series
ser = series(WorldClim{Climate}, :tavg; month=1:12)

# Run for 12 months of 2017, on a CUDA GPU
@time rates = mapgrowth((host=host_model, para=para_model);
    series=ser,
    tspan=1:12,
)

# Save tifs for each layer
for t in 1:12
    path = joinpath("output/growthrates/")
    mkpath(joinpath(path, "host"))
    mkpath(joinpath(path, "parasite"))
    tpad = lpad(t, 2, '0')
    write(joinpath(path, "host/$tpad.tif"), GDALarray, rates[:host][Ti(t)])
    write(joinpath(path, "parasite/$tpad.tif"), GDALarray, rates[:para][Ti(t)])
end

#### Simple Growth Model #####################

using GrowthMaps, Unitful, GeoData, HDF5, Dates, CUDA, Random, KernelAbstractions
using GeoData: Between
CUDA.allowscalar(false)

# Weird: Define julia object to hold parameters to for the simple growth model 
# a struct with parameters, that inherits from GrowthModel *** more about this
struct IntrinsicGrowth{Ic,U,If,D} <: GrowthModel
    intercept::Ic
    up::U
    inflection::If
    down::D
end
# And a Growthmaps.rate method to calculate growth for a specific temperature
# Weird: we are adding a method to a function in the GrowthMaps mackage
# for the struct we have defined above.
GrowthMaps.rate(m::IntrinsicGrowth, temp) = 
    gr(m.intercept, m.up, m.inflection, m.down, temp)
gr(a, b, c, d, x) = x < c ? a + b*x  : a + b*c - d*(x - c)  
# Parametrise the struct
host_growth = IntrinsicGrowth(-14.6, 0.05u"K^-1", 30u"°C" |> u"K", 0.1u"K^-1")
para_growth = IntrinsicGrowth(-14.6, 0.05u"K^-1", 26u"°C" |> u"K", 0.1u"K^-1")
# Specify the SMAP layer `surface_temp` with units Kelvin
# TODO: use simplified data in the repo
host_model = Model(Layer(:surface_temp, u"K", host_growth))
para_model = Model(Layer(:surface_temp, u"K", para_growth))

# Plot this with scatter, if there is room
temps = collect(-15.0u"°C":0.1u"K":40.0u"°C")

# Get time series of SMAP files to load
# days = 1:31 # all
days = 1 # one day a month
series = SMAPseries("/home/raf/Data/SMAP/SMAP_L4_SM_gph_v4")[Where(t -> dayofmonth(t) in days)]

# Run for 12 months of 2017, on a CUDA GPU
@time rates = mapgrowth(host_model, para_model;
    series=series,
    tspan=DateTime(2017, 1):Month(1):DateTime(2017, 12),
    arraytype=CuArray, # Run on a CUDA GPU
)

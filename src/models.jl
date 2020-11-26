#### Simple Growth Model #####################

using GrowthMaps, CUDA, Plots, Unitful, GeoData, HDF5, Dates, DynamicGridsInteract, DynamicGrids
using DynamicGrids: auxval
const DG = DynamicGrids
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
growthmodel = IntrinsicGrowth(-14.6, 0.05u"K^-1", 30u"°C" |> u"K", 0.1u"K^-1")
# Specify the SMAP layer `surface_temp` with units Kelvin
# TODO: use simplified data in the repo
model = Model(Layer(:surface_temp, u"K", growthmodel))

# Plot this with scatter, if there is room
temps = collect(-15.0u"°C":0.1u"K":40.0u"°C")

# Get time series of SMAP files to load
# days = 1:31 # all
days = 1 # one day a month
series = SMAPseries("/home/raf/Data/SMAP/SMAP_L4_SM_gph_v4")[Where(t -> dayofmonth(t) in days)]

# Run for 12 months of 2017, on a CUDA GPU
@time growthrates = mapgrowth(model;
    series=series,
    tspan=DateTime(2017, 1):Month(1):DateTime(2017, 12),
    arraytype=CuArray, # Run on a CUDA GPU
)

aus = Lon(Between(113.3402, 153.9523)), Lat(Between(-43.62234, -10.65125))
growthrates1 = growthrates[aus...] |> 
               # gr -> set(gr, Lat(Base.OneTo(length(dims(growthrates, Lat))), NoIndex())) |>
               # gr -> set(gr, Lon(Base.OneTo(length(dims(growthrates, Lon))), NoIndex())) |>
               # gr -> set(gr, NoName()) |> 
               gr -> replace_missing(gr, 0.0f0) |> 
               gr -> permutedims(gr, (Lat, Lon, Ti))

# ------------------
# Implications, useful for all of SDM modelling as a side effect

# Growth model for the dispersal simulation,
# using the growth rates we have just generated:
# Here we use a pre-defined growth model LogisticGrowthMap from Dispersal.jl
using Dispersal, ColorSchemes, BenchmarkTools
import Dispersal: applyrule, aux
carrycap = 1f9
growth = LogisticGrowthMap{:H,:H}(
    auxkey=Val(:rH),
    carrycap=carrycap,
    timestep=Day(1),
)

#### Host dispersal model #####################
λ = Param(0.125f0; bounds=(0.01, 0.3))
radius = 1
hood = DispersalKernel{radius}(; formulation=ExponentialKernel(λ), cellsize=1.0f0)
localdisp = InwardsPopulationDispersal{:H,:H}(hood)

# Output 
hostpop = zero(growthrates1[Ti(1)])
sydney = Lat(Near(-33)), Lon(Between(150.0, 151.0))
cairns = Lat(Between(-16, -19)), Lon(Between(134.0, 147.0))
hostpop[cairns...] .= carrycap

tspan = DateTime(2017, 1):Day(7):DateTime(2030, 1)
mask = boolmask(growthrates1[Ti(1)])


#=
This time we define our own model from scratch, 
using a shorthand model definition style - which acheives brevity 
with some loss of reusability. Well run this on a CPU,
as generating random numbers on a GPU is one of the things
that doesn't "just work" on GPU, and needs a different approach
to the simple method we use here.
=#


### Wind dispersal model #####################################
wind = Manual{:H,:H}() do data, index, state
    # Ignore empty cells
    state > zero(state) || return nothing
    # Randomise a dispersal event
    rand() < 0.01 || return nothing
    # Randomise a destination
    spotrange = 20
    rnge = -spotrange:spotrange
    jump = (rand(rnge), rand(rnge))
    # Make sure the destination is on the grid
    jumpdest, is_inbounds = inbounds(jump .+ index, gridsize(data), RemoveOverflow())
    # Update spotted cell if it's on the grid
    if is_inbounds 
        @inbounds add!(data[:H], state / 5, jumpdest...)
        @inbounds sub!(data[:H], state / 5, index...)
    end
    return nothing
end

using DynamicGridsGtk
viewer = AuxCopy{:rH}(Val(:rH))
output = GtkOutput((; H=hostpop, rH=zero(growthrates1[Ti(1)])); 
    filename="sim.gif", 
    aux=(; rH=growthrates1),
    mask=parent(mask),
    tspan=tspan, 
    fps=1000,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, maximum(growthrates1)),
    scheme=(ColorSchemes.autumn, ColorSchemes.inferno),
)
sim!(output, (localdisp, growth, viewer))
sim!(output, (wind, growth, viewer))
sim!(output, (wind, localdisp, growth, viewer))


#### Host/Parasite growth model #####################
struct HostParasiteGrowth{R,W,M,P,Hs,Hh} <: CellRule{R,W}
    γ::M
    ρ::P
    Hsat::Hs
    Hhalf::Hh
end

hollingII(H, β, h) = β * H / (1 + h * β * H)
encounter_rate(Hsat, Hhalf, h) = (Hsat - 2Hhalf) / (h * Hsat * Hhalf)

function DynamicGrids.applyrule(data, rule::HostParasiteGrowth, (H, P), index)
    # Parameters
    γ = rule.γ; ρ = rule.ρ; Hsat = rule.Hsat; Hhalf = rule.Hhalf
    # Growth rates for this cell
    rH = auxval(data, :rH, index...)
    rP = auxval(data, :rP, index...)
    # Equation
    h = ρ * rP
    β = encounter_rate(Hsat, Hhalf, h)
    R = hollingII(H, β, h)
    δ = (rP + γ) / R
    dHdt = rH * H - R * P
    dPdt = R * δ * P - γ * P
    H + dHdt, P + dPdt
end

γ = Param(0.024; bounds=(0.01, 0.5)) 
ρ = Param(63.8; bounds=(0.0, 100.0)) 
Hsat = Param(1e7; bounds=(0.0, 1e10)) 
Hhalf = Param(1e6; bounds=(0.0, 1e10))
host_para_growth = HostParasiteGrowth{Tuple{:H,:P},Tuple{:H,:P}}(γ, ρ, Hsat, Hhalf)
λ = 0.05
radius = 1
para_hood = DispersalKernel{radius}(; formulation=ExponentialKernel(λ))
para_localdisp = InwardsPopulationDispersal{:P,:P}(para_hood)
parapop = (a -> (r = rand(Float32); r < 0.0001 ? carrycap : 0.0f0)).(hostpop)
ruleset = Ruleset(wind, host_para_growth)
output = GtkOutput((H=hostpop, P=parapop); 
    filename="sim.gif", 
    ruleset=ruleset,
    aux=(; rH=growthrates1, rP=growthrates1),
    mask=parent(mask),
    tspan=tspan, 
    fps=10,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, carrycap),
    scheme=(ColorSchemes.autumn, ColorSchemes.autumn),
)
display(output)
sim!(output, ruleset)

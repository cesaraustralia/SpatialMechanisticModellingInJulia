using GeoData, ArchGDAL, Dates

path = "/home/raf/.julia/dev/Leafminer/output/growthrates/huidobrensis"
gr_slices = map(readdir(path; join=true)) do path
    GDALarray(path; mappedcrs=EPSG(4326))[Band(1)] |>
        a->permutedims(a, (Lat, Lon))
end
growthtimespan = DateTime(2017,1):Month(1):DateTime(2017,12)
gr = cat(gr_slices...; dims=Ti(growthtimespan))
aus = Lon(Between(113.3402, 153.9523)), Lat(Between(-43.62234, -10.65125))
rH = gr[aus...] |> gr -> replace_missing(gr, 0.0f0) |> gr -> permutedims(gr, (Lat, Lon, Ti))
rP = rH


### Growth model #####################
using Dispersal, ColorSchemes, BenchmarkTools, DynamicGridsInteract, DynamicGrids
const DG = DynamicGrids

carrycap = 1f9
growth = LogisticGrowth{:H}(
    rate=Aux(:rH),
    carrycap=carrycap,#Param(carrycap, bounds=(0.0, 10carrycap)),
    timestep=Day(1),
)

# Output 
hostpop = zero(rH[Ti(1)])
cairns = Lat(Between(-17, -20)), Lon(Between(144.0, 145.0))
hostpop[cairns...] .= carrycap
tspan = DateTime(2020, 1):Day(7):DateTime(2044, 1)
mask = boolmask(rH[Ti(1)])

alleeH = AlleeExtinction{:H}(; minfounders=20.0)

# output = ElectronOutput((H=hostpop, rH=zero(hostpop)); 
output = GtkOutput((H=hostpop, rH=zero(hostpop)); 
    ruleset = Ruleset(growth, viewer),
    filename="sim.gif", 
    aux=(; rH=rH),
    mask=parent(mask),
    tspan=tspan, 
    store=true,
    fps=50,
    minval=(0.0f0, 0.0), 
    maxval=(carrycap, maximum(rH)),
    scheme=(ColorSchemes.autumn, ColorSchemes.inferno),
)
sim!(output, (growth, viewer)) 

output = ArrayOutput((; H=fill(carrycap, size(hostpop))); 
    aux=(; rH=rH),
    mask=mask,
    tspan=tspan, 
)
sim!(output, growth) 


#### Host dispersal model #####################
λ = Param(0.25f0; bounds=(0.01, 0.3))
r = 1
hood = DispersalKernel{r}(; formulation=ExponentialKernel(λ), cellsize=1.0f0)
localdisp = InwardsDispersal{:H}(hood)


### Wind dispersal model #####################################
wind = SetCell{:H}() do data, I, state 
    state > zero(state) || return nothing # Ignore empty cells
    rand() < 0.01 || return nothing # Randomise a dispersal event
    jump = (rand(-20:20), rand(-20:20)) # Randomise a destination
    # Make sure the destination is on the grid
    jumpdest, is_inbounds = inbounds(jump .+ I, data)
    # Update spotted cell if it's on the grid
    if is_inbounds 
        @inbounds add!(data[:H], state / 10, jumpdest...)
        @inbounds sub!(data[:H], state / 10, I...)
    end
    return nothing
end

# GPU wind model
# We can't use `rand` directly inside a GPU kernel, as of Dec 2020.
# To get around this shortcoming we add a grid to hold random numbers and 
# generate them from outside a kernel using a `SetGrid` model that calls a function
# on the entire grid. All the work is still done on the GPU.
gpu_wind = SetNeighbors{Tuple{:rand,:H}}() do data, hood, I, rand, H
    # Ignore empty cells
    H > zero(H) || return nothing
    # Randomise a dispersal event using cell value of rand,
    # calculated in the rand Grid rule
    rand < 0.02f0 || return nothing
    # Randomise a destination
    x = 10oneunit(H)
    # Get random values for jump from the surrounding grid
    randi = (I[1], I[2] + 1)
    randj = (I[1] + 1, I[2] + 1)
    jump = unsafe_trunc(Int, data[:rand][randi...] * 2x - x), 
           unsafe_trunc(Int, data[:rand][randj...] * 2x - x)
    # Make sure the destination is on the grid
    jumpdest, is_inbounds = inbounds(jump .+ I, data[:H])
    # Update spotted cell if it's on the grid
    if is_inbounds 
        @inbounds add!(data[:H], H / 10, jumpdest...)
        @inbounds sub!(data[:H], H / 10, I...)
    end
    return nothing
end

_rand!(w::CuArray) = CUDA.rand!(w)
_rand!(w) = rand!(w)

randomgrid = SetGrid{Tuple{},:rand}() do w
    _rand!(parent(w))
end

using DynamicGridsGtk
output = GtkOutput((; H=hostpop); 
    ruleset = Ruleset(wind, localdisp, growth),
    filename="sim.gif", 
    aux=(; rH=rH),
    mask=parent(mask),
    tspan=tspan, 
    store=true,
    fps=50,
    minval=(0.0f0, 0.0f0, 0.0f0), 
    maxval=(carrycap, maximum(rH), 1.0f0),
    scheme=(ColorSchemes.autumn, ColorSchemes.inferno, ColorSchemes.inferno),
)


sim!(output, (wind, growth));

# DynamicGrids.savegif("localdisp.gif", output)
# sim!(output, (wind, growth, viewer))
# DynamicGrids.savegif("winddisp.gif", output)
# sim!(output, (wind, localdisp, growth, viewer))
# DynamicGrids.savegif("combinded.gif", output)


#### Host/Parasite growth model #####################
struct HostParasiteGrowth{R,W,M,P,Hs,Hh} <: CellRule{R,W}
    γ::M
    ρₕ::P
    Hsat::Hs
    Hhalf::Hh
end

typeII(H, β, h) = β * H / (1 + h * β * H)
# encounter_rate(Hsat, Hhalf, h) = (Hsat - 2Hhalf) / (h * Hsat * Hhalf)
encounter_rate(Hhalf, h) = 1 / (h * Hhalf)

@inline function DynamicGrids.applyrule(data, rule::HostParasiteGrowth, (H, P), I)
    # Parameters
    γ = rule.γ; Hhalf = rule.Hhalf
    # Growth rates for this cell
    rH = get(data, Aux(:rH), I...)
    rP = get(data, Aux(:rP), I...)
    # Equation
    h = 1 / (2.161057 * exp(9.050459 * rP))
    β = encounter_rate(Hhalf, h)
    R = typeII(H, β, h)
    δ = h * (rP - γ)
    dHdt = rH * H - R * P
    dPdt = δ * R * P * rP - γ * P
    # if H > zero(H)
        # @show H P rH rP h β R δ dHdt dPdt
        # println()
    # end
    H + dHdt, P + dPdt
end

sd = DynamicGrids.SimData(output.extent, Ruleset());
sd = DynamicGrids._updatetime(sd, 1);
DynamicGrids.applyrule(sd, host_para_growth, (2000000.0, 20000.0), (96, 176))

γ = Param(0.024; bounds=(0.01, 0.5)) 
ρ = Param(63.8; bounds=(0.0, 100.0)) 
Hsat = Param(1e7; bounds=(0.0, 1e10)) 
Hhalf = Param(1e2; bounds=(0.0, 1e10))
host_para_growth = HostParasiteGrowth{Tuple{:H,:P},Tuple{:H,:P}}(γ, ρ, Hsat, Hhalf)

λ = Param(0.25f0, bounds=(0.0, 0.3))
r = 1
para_hood = DispersalKernel{r}(; formulation=ExponentialKernel(λ))
para_localdisp = InwardsDispersal{:P}(para_hood)
alleeP = AlleeExtinction{:P}(; minfounders=20.0)

parapop = (a -> (r = rand(Float32); r < 0.01 ? carrycap : 0.0f0)).(hostpop)
# # fill!(parapop, 1e8)
# ruleset = Ruleset(gpu_wind, randomgrid, localdisp, para_localdisp, host_para_growth; proc=CuGPU())
# ruleset = Ruleset(wind, localdisp, para_localdisp, host_para_growth)
# # output = ElectronOutput((H=hostpop, P=parapop, rand=zero(hostpop)); 
# output = GtkOutput((H=hostpop, P=parapop, rand=zero(hostpop)); 
#     filename="sim.gif", 
#     ruleset=ruleset,
#     aux=(; rH=rH, rP=rP),
#     mask=parent(mask),
#     tspan=tspan, 
#     fps=10,
#     minval=(0.0f0, 0.0f0, 0.0f0), 
#     maxval=(carrycap, carrycap, 1.0f0),
#     scheme=(ColorSchemes.autumn, ColorSchemes.autumn, ColorSchemes.autumn),
# )
# display(output)
# sim!(output, ruleset)


copy_rH = CopyTo{:rH}(Aux{:rH}())
# output = ElectronOutput((H=hostpop, P=parapop, rand=zero(hostpop)); 
using DynamicGridsGtk
output = GtkOutput((H=hostpop, P=parapop, rH=zero(parapop)); 
    filename="sim.gif", 
    aux=(; rH=rH, rP=rP),
    mask=parent(mask),
    tspan=tspan,
    fps=100,
    store=true,
    minval=(0.0f0, 0.0f0, 0.f0), 
    maxval=(carrycap, carrycap, maximum(rH)),
    scheme=(ColorSchemes.autumn, ColorSchemes.autumn, ColorSchemes.inferno),
)

sim!(output, wind, alleeH, alleeP, localdisp, para_localdisp, host_para_growth, copy_rH)

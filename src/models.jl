
using Dispersal, Dates

### Growth rule #######################################
carrycap = 1f9
growth = LogisticGrowth{:H}(
    rate=Aux(:rH),
    carrycap=Param(carrycap, bounds=(1f0, 10carrycap)),
    timestep=Day(1),
)

### Allee extinction rule #############################
allee = AlleeExtinction{:H}(; minfounders=Param(20.0f0; bounds=(1f0, 1e3)))

#### Host local dispersal rule ########################
λ = Param(0.25f0; bounds=(0.01f0, 0.3f0))
r = 1
hood = DispersalKernel{r}(; formulation=ExponentialKernel(λ), cellsize=1.0f0)
localdisp = InwardsDispersal{:H}(hood)


### Wind dispersal rule ###############################
wind = SetCell{:H}() do data, I, state 
    state > zero(state) || return nothing # Ignore empty cells
    rand() < 0.01f0 || return nothing # Randomise a dispersal event
    jump = rand(-20:20), rand(-20:20) # Randomise a destination
    # Make sure the destination is on the grid
    jumpdest, is_inbounds = inbounds(jump .+ I, data)
    # Update spotted cell if it's on the grid
    if is_inbounds 
        @inbounds add!(data[:H], state / 10, jumpdest...)
        @inbounds sub!(data[:H], state / 10, I...)
    end
    return nothing
end

# GPU wind rule
# We can't use `rand` directly inside a GPU kernel, as of Dec 2020.
# To get around this shortcoming we add a grid to hold random numbers and 
# generate them from outside a kernel using a `SetGrid` rule that calls a function
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

# _rand!(w::CuArray) = CUDA.rand!(w)
_rand!(w) = rand!(w)

randomgrid = SetGrid{Tuple{},:rand}() do w
    _rand!(parent(w))
end


#### Host/Parasite growth rule #####################
struct HostParasite{R,W,G,Hh,A,B} <: CellRule{R,W}
    γ::G
    Hhalf::Hh
    a::A
    b::B
end

typeII(H, β, h) = β * H / (1 + h * β * H)
encounter_rate(Hhalf, h) = 1 / (h * Hhalf)
# encounter_rate(Hsat, Hhalf, h) = (Hsat - 2Hhalf) / (h * Hsat * Hhalf)

@inline function Dispersal.applyrule(data, rule::HostParasite, (H, P), index)
    # Parameters
    γ = rule.γ; Hhalf = rule.Hhalf; a = rule.a; b = rule.b
    # Growth rates for this cell
    rH = get(data, Aux(:rH), index...)
    rP = get(data, Aux(:rP), index...)
    # Equation
    h = 1 / (a * exp(b * rP))
    β = encounter_rate(Hhalf, h)
    R = typeII(H, β, h)
    δ = h * (rP - γ)
    dHdt = rH * H - R * P
    dPdt = δ * R * P * rP - γ * P
    H + dHdt, P + dPdt
end

# Define the object
# ρ = Param(63.8; bounds=(0.0, 100.0)) 
# Hsat = Param(1e7; bounds=(0.0, 1e10)) 
γ = Param(0.024f0; bounds=(0.01f0, 0.5f0)) 
Hhalf = Param(1f9; bounds=(0.0f0, 1f10))
a = Param(2.161057f0; bounds=(1f-7, 1f1))
b = Param(9.050459f0; bounds=(1f-7, 1f2))
host_parasite = HostParasite{Tuple{:H,:P}}(γ, Hhalf, a, b)

#### Parasite local dispersal rule ##################
λ = Param(0.25f0, bounds=(0.0f0, 0.3f0))
r = 1
localdispP = InwardsDispersal{:P}(
    neighborhood=DispersalKernel{r}(; formulation=ExponentialKernel(λ))
)
alleeP = AlleeExtinction{:P}(; minfounders=Param(20.0f0; bounds=(1f0, 1e3)))

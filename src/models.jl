using Dispersal, Dates

### Growth rule #######################################
carrycap = 1f9
growth = LogisticGrowth{:H}(
    rate=Aux(:rH), carrycap=carrycap, timestep=Day(1), nsteps_type=Float32,
)

### Allee extinction rule #############################
allee = AlleeExtinction{:H}(; minfounders=20.0f0)

#### Host local dispersal rule ########################
λ = 0.25f0
radius = 1
hood = DispersalKernel{radius}(; formulation=ExponentialKernel(λ), cellsize=1.0f0)
localdisp = InwardsDispersal{:H}(hood)

### Wind dispersal rule ###############################
wind = SetCell{:H}() do data, I, N 
    N > zero(N) || return nothing # Ignore empty cells
    rand() < 0.01f0 || return nothing # Randomise a dispersal event
    dest = I .+ (rand(-20:20), rand(-20:20)) # Randomise a destination
    # Update spotted cell, if it's on the grid
    if isinboknds(dest, data)
        @inbounds add!(data[:H], N / 10, dest...)
        @inbounds sub!(data[:H], N / 10, I...)
    end
    return nothing
end

# GPU wind rule
# We can't use `rand` directly inside a GPU kernel, as of February 2021.
# To get around this shortcoming we add a grid to hold random numbers and 
# generate them from outside a kernel using a `SetGrid` rule that calls a function
# on the entire grid. All the work is still done on the GPU.
gpu_wind = SetNeighbors{Tuple{:rand,:H}}() do data, hood, I, rand, H
    # Ignore empty cells
    H > zero(H) || return nothing
    # Randomise a dispersal event using cell value of rand
    rand < 0.02f0 || return nothing
    x = 10oneunit(H) # Randomise a destination
    # Get random values for jump from the surrounding grid
    randi, randj = (I[1], I[2] + 1), (I[1] + 1, I[2] + 1)
    dest = I .+ (unsafe_trunc(Int, data[:rand][randi...] * 2x - x), 
                 unsafe_trunc(Int, data[:rand][randj...] * 2x - x))
    # Update spotted cell if it's on the grid
    if isinbounds(dest, data[:H])
        @inbounds add!(data[:H], H / 10, dest...)
        @inbounds sub!(data[:H], H / 10, I...)
    end
    return nothing
end

_rand!(A) = rand!(A)

randomgrid = SetGrid{Tuple{},:rand}() do g
    _rand!(parent(g))
end

#### Host/Parasite growth rule #####################
struct Parasitism{R,W,RP,G,Hh,A,B} <: CellRule{R,W}
    rP::RP
    γ::G
    Hhalf::Hh
    a::A
    b::B
end

function Dispersal.applyrule(data, rule::Parasitism, (H, P), index)
    # Parameters
    γ = rule.γ; Hhalf = rule.Hhalf; a = rule.a; b = rule.b
    # Parasite maximum growth rate for this cell
    rP = get(data, rule.rP, index...) * 7
    # Equation
    h = 1 / (a * exp(b * rP))
    β = 1 / (h * Hhalf) # Encounter rate
    R = β * H / (1 + h * β * H) # Holling type II response
    δ = h * (rP - γ)
    dHdt = -R * P
    dPdt = δ * R * P * rP - γ * P
    return H + dHdt, P + dPdt
end

# Construct a Parasitism rule
parasitism = Parasitism{Tuple{:H,:P}}(Aux{:rP}(), 0.2f0, 1f8, 2.0f0, 5.0f0)

#### Parasite local dispersal rule ##################
λ = 0.15f0
radius = 1
localdispP = InwardsDispersal{:P}(
    neighborhood=DispersalKernel{radius}(; formulation=ExponentialKernel(λ))
)
alleeP = AlleeExtinction{:P}(; minfounders=20.0f0)

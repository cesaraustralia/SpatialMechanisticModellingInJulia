using Dispersal, Dates

### Growth rule #######################################
carrycap = 1f9
growth = LogisticGrowth{:H}(
    rate=Aux(:rH), carrycap=carrycap, timestep=Day(1), nsteps_type=Float32,
)

### Allee extinction rule #############################
allee = AlleeExtinction{:H}(; minfounders=20.0f0)

#### Host local dispersal rule ########################
hood = DispersalKernel(; 
    radius=1, formulation=ExponentialKernel(0.25f0), cellsize=1.0f0
)
localdisp = InwardsDispersal{:H}(hood)

### Wind dispersal rule ###############################
wind = SetCell{:H}() do data, I, N 
    N > zero(N) || return nothing # Ignore empty cells
    rand() < 0.01f0 || return nothing # Randomise a dispersal event
    dest = I .+ (rand(-20:20), rand(-20:20)) # Randomise a destination
    # Update spotted cell, if it's on the grid
    if isinbounds(dest, data)
        @inbounds add!(data[:H], N / 10, dest...)
        @inbounds sub!(data[:H], N / 10, I...)
    end
    return nothing
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
λ = 0.18f0
radius = 1
localdispP = InwardsDispersal{:P}(
    neighborhood=DispersalKernel{radius}(; formulation=ExponentialKernel(λ))
)
alleeP = AlleeExtinction{:P}(; minfounders=20.0f0)

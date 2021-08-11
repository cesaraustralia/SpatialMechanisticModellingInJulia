using Dispersal, Dates

#### Parasitoid growth rule object and method #####################
struct Parasitism{R,W,GR,M,HS,A,B} <: CellRule{R,W}
    rate::GR
    mortality::M
    halfsat::HS
    a::A
    b::B
end

function Dispersal.applyrule(data, rule::Parasitism, (H, P), I)
    # Parameters
    γ = rule.mortality; halfsat = rule.halfsat; a = rule.a; b = rule.b
    # Parasite maximum growth rate for this cell, may vary spatially
    rP = get(data, rule.rate, I...) * 7
    # Equation
    h = 1 / (a * exp(b * rP))
    β = 1 / (h * halfsat) # Encounter rate
    R = β * H / (1 + h * β * H) # Holling type II response
    δ = h * (rP - γ)
    dHdt = -R * P
    dPdt = δ * R * P * rP - γ * P
    return H + dHdt, P + dPdt
end

function define_rules()

    ### Growth rule #######################################
    carrycap = 1f9
    growth = LogisticGrowth{:H}(
        rate=Aux(:rH), carrycap=carrycap, timestep=Day(1), 
        nsteps_type=Float32, # for type stability and performance with a Float32 grid
    )

    ### Allee extinction rule #############################
    allee = AlleeExtinction{:H}(; minfounders=20.0f0)

    #### Host local dispersal rule ########################
    localdisp = InwardsDispersal{:H}(;
        radius=1, formulation=ExponentialKernel(0.25f0), cellsize=1.0f0,
    )

    ### Wind dispersal rule ###############################
    wind = SetCell{:H}() do data, H, I 
        H > zero(H) || return nothing # Ignore empty cells
        rand() < 0.01f0 || return nothing # Randomise a dispersal event
        # Randomise a destination
        dest = I .+ (rand(-20:20), rand(-20:20))
        if isinbounds(data, dest)
            # Update dest cell and current cell, if dest is on the grid
            @inbounds add!(data[:H], H / 10, dest...)
            @inbounds sub!(data[:H], H / 10, I...)
        end
        return nothing
    end

    # Construct a Parasitism rule
    parasitism = Parasitism{Tuple{:H,:P}}(Aux(:rP), 0.2f0, 1f8, 2.0f0, 5.0f0)

    #### Parasitoid local dispersal rules ##################
    localdisp_p = InwardsDispersal{:P}(; 
        radius=1, formulation=ExponentialKernel(; λ=0.18f0), cellsize=1.0f0
    )
    allee_p = AlleeExtinction{:P}(; minfounders=20.0f0)

    return (; growth, allee, localdisp, wind, parasitism, localdisp_p, allee_p, carrycap)
end


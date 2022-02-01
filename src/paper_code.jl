using Dispersal, Dates
growth = LogisticGrowth{:H}(
    rate=Aux(:rH), carrycap=1e9, timestep=Day(1)
)


localdisp = InwardsDispersal{:H}(
    radius=1, formulation=ExponentialKernel(λ = 0.25), cellsize=1.0
)
allee = AlleeExtinction{:H}(; minfounders=20.0)


wind = SetCell{:H}() do data, H, I
    H > zero(H) || return nothing # Ignore empty cells
    rand() < 0.01 || return nothing # Randomise a rare dispersal event
    dest = I .+ (rand(-20:20), rand(-20:20)) # Randomise a destination
    if isinbounds(data, dest) # Check i# f dest is on the grid
        add!(data[:H], H / 10, dest...) # Add population to dest cell
        sub!(data[:H], H / 10, I...) # Subtract from current cell
    end
    return nothing
end

# To generate the data variables see Appendix S3
output = ArrayOutput((H=init_h,);
    tspan=DateTime(2020, 1):Week(1):DateTime(2022, 1),
    aux=(rH=rate_h,), mask=mask,
) 
# Run the rules for different model combinations 
sim!(output, (localdisp, allee, growth))       # local
sim!(output, (wind, allee, growth))            # wind
sim!(output, (wind, localdisp, allee, growth)) # combined
using Plots
plot(output[end][:H]; color=:inferno) # Plot the last frame of a simulation


struct Parasitism{R,W,GR,M,HS,A,B} <: CellRule{R,W}
    rate::GR
    mortality::M
    halfsat::HS
    a::A
    b::B
end

function Dispersal.applyrule(data, rule::Parasitism, (H, P), I)
    γ = rule.mortality; halfsat = rule.halfsat; a = rule.a; b = rule.b
    rP = get(data, rule.rate, I) * 7 # Growthrate for the cell/timestep
    h = 1 / (a * exp(b * rP)) # Handling time
    β = 1 / (h * halfsat) # Encounter rate
    R = β * H / (1 + h * β * H) # Holling type II response
    δ = h * (rP - γ) # Conversion efficiency
    dHdt = -R * P # Change in host population
    dPdt = δ * R * P * rP - γ * P # Change in parasitoid population
    return (H + dHdt, P + dPdt) # Updated populations for the cell
end

parasitism = Parasitism{Tuple{:H,:P}}(Aux(:rP), 0.2, 1e8, 2.0, 5.0)

localdisp_p = InwardsDispersal{:P}(; 
    radius=1, formulation=ExponentialKernel(; λ=0.18), cellsize=1.0
)
allee_p = AlleeExtinction{:P}(; minfounders=20.0)
output = GifOutput((H=init_h, P=init_p); 
    filename="host_parasitoid.gif", 
    tspan=DateTime(2020, 1):Week(1):DateTime(2030, 1),
    aux=(rH=rate_h, rP=rate_p,),
    mask=mask,
    minval=(0.0, 0.0),
    maxval=(1e9, 1e9),
) 
# Run and save as "host_parasitoid.gif"
sim!(output, (wind, localdisp, allee, growth, localdisp_p, allee_p, parasitism))

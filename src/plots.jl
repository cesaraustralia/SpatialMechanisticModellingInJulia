
### Plots ######################################################################

using Shapefile, Plots, Colors
using Plots: px

# Plot setup
shp = Shapefile.Handle("/home/raf/.julia/dev/MethodsPaper/data/ausborder_polyline.shp")
shape_opts = (lw=0.5, color=RGB(0.3))
plot_opts = (; 
    margin=0px, 
    foreground_color=RGB(0.3),
    tickfontsize=7, 
    title="", 
    xlims=(112, 155), ylims=(-44, -8)
)

# Growth and environment plots
growthplot_opts = (; plot_opts..., tickfontcolor=RGB(0.5))
p1 = plot(series[Ti(1)][:surface_temp][aus...]; 
    colorbar_title="Surface temperature", c=:acton, xtickfontcolor=RGB(1.0), growthplot_opts...
)
p2 = plot(growthrates[Ti(1), aus...]; 
    colorbar_title="Intrinsic growth", c=:thermal, growthplot_opts...
)
map(p -> plot!(p, shp; shape_opts...), (p1, p2))
growthplots = plot(p1, p2; layout=(2, 1), link=:x, size=(600, 900))
# savefig("output/growthplots.png")

# Simulation plots
init = (; H=hostpop)
kwargs = (; 
    aux=(; rH=rH, rP=rP), 
    mask=parent(mask),
    tspan=tspan
)
wind_output = ArrayOutput(init; kwargs...)
localdisp_output = ArrayOutput(init; kwargs...)
combined_output = ArrayOutput(init; kwargs...)
growth_output = ArrayOutput((; H=copy(hostpop) .= 1e7); kwargs...)

sim!(wind_output, wind, allee, growth)
sim!(localdisp_output, localdisp, allee, growth)
sim!(combined_output, localdisp, wind, allee, growth)
sim!(growth_output, allee, growth)

simplot_opts = (; 
    xtickfontcolor=RGB(0.5), 
    clims=(0, carrycap),
    showaxis=false,
    legend=:none, plot_opts...
)
plot_ti = Ti(Near(DateTime(2023, 1)))
p3 = plot(localdisp_output[plot_ti][:H]; xguide="Local", simplot_opts..., ytickfontcolor=RGB(0.5))
p4 = plot(wind_output[plot_ti][:H]; xguide="Wind", simplot_opts..., ytickfontcolor=RGB(1.0))
p5 = plot(combined_output[plot_ti][:H]; xguide="Combined", simplot_opts..., ytickfontcolor=RGB(1.0))
map(p -> plot!(p, shp; shape_opts...), (p3, p4, p5))
simplots = plot(p3, p4, p5; layout=(1, 3), size=(1000, 400))
savefig("output/simplots.png")

# Location line plot
locations = (
#    Adelaide=(-34.901608, 138.601547),
    Brisbane=(-27.436190, 152.990588),
    Cairns=(-16.937281, 145.747709),
    Darwin=(-12.4634, 130.8456),
    Hobart=(-42.881742, 147.323879),
    Melbourne=(-37.805896, 144.959527),
    Perth=(-31.9505, 115.8605),
    Sydney=(-33.839943, 151.006101),
)
plotspan = Between(DateTime(2020), DateTime(2022))
loc_populations = map(locations) do loc
    coords = Lat(Contains(loc[1])), Lon(Contains(loc[2]))
    map(f -> f[:H][coords...] + 1, growth_output[plotspan])
end
map(f -> f[:H][100, 100], growth_output)
growth_output[end][:H][100, 100]
plot(growth_output[end][:H])
theme(:vibrant)

growthplot = plot(Date.(index(growth_output[plotspan], Ti)), reduce(hcat, loc_populations); 
    labels=reshape([string.(keys(locations))...], (1, length(locations))), 
    lalpha=0.7,
    xguide="Time", 
    yguide="Population",
    yscale=:log10,
    ylims=(2.0, 1e9),
    legend=:bottomright,
    margin=20px,
    size=(400, 400),
)

savefig("output/growthplot.png")



init = (H=hostpop, P=parapop)
host_para_output = ArrayOutput(init; kwargs...)
sim!(host_para_output, localdisp, wind, host_para_growth)

p7 = plot(host_para_output[300][:H]; yguide="Host-Parasite", simplot_opts..., xtickfontcolor=RGB(1.0))
p8 = plot(host_para_output[600][:H]; yguide="Host-Parasite", simplot_opts..., xtickfontcolor=RGB(1.0))
hp_plots = plot(p7, p8; layout=(2, 1), size=(600, 600))

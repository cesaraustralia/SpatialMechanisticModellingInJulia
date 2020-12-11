
### Plots ######################################################################

using Shapefile, Plots, Colors
using Plots: px

# Plot setup
shp = Shapefile.Handle("/home/raf/.julia/dev/MethodsPaper/data/ausborder_polyline.shp")
shape_opts = (lw=1.0, color=RGB(0.3))
plot_opts = (; 
    bottom_margin=0px, 
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
init = (; H=copy(hostpop) .= 10000)
growth_output = ArrayOutput(; kwargs...)

sim!(wind_output, wind, growth)
sim!(localdisp_output, localdisp, growth)
sim!(combined_output, localdisp, wind, growth)
sim!(growth_output, growth)

output = ArrayOutput(init; 
    aux=(; rH=hH, rP=rH), 
    mask=boolmask(,
    tspan=DateTime(2020, 1):Day(7):DateTime(2024, 1)
)

simplot_opts = (; 
    yguide_position=:right,
    ytickfontcolor=RGB(0.5), 
    clims=(0, carrycap),
    legend=:none, plot_opts...
)
plotti = Ti(DateTime(2025, 1))
p3 = plot(localdisp_output[plotti][:H]; yguide="Local", simplot_opts..., xtickfontcolor=RGB(1.0))
p4 = plot(wind_output[plotti][:H]; yguide="Wind", simplot_opts..., xtickfontcolor=RGB(1.0))
p5 = plot(combined_output[plotti][:H]; yguide="Combined", simplot_opts..., xtickfontcolor=RGB(0.5))
map(p -> plot!(p, shp; shape_opts...), (p3, p4, p5))
simplots = plot(p3, p4, p5; layout=(3, 1), size=(600, 900))
# savefig("output/simplots.png")

locations = (
    Melbourne=(-37.805896, 144.959527),
    Mildura=(-34.219504, 142.130864),
    Coffs_Harbour=(-30.287245, 153.092991),
    Sydney=(-33.839943, 151.006101),
    Adelaide=(-34.901608, 138.601547),
    Port_Augusta=(-32.466201, 137.813850),
    Devonport=(-41.180545, 146.314887),
    Hobart=(-42.881742, 147.323879),
    Brisbane=(-27.436190, 152.990588),
    Cairns=(-16.937281, 145.747709),
    Perth=(-31.9505, 115.8605),
    Geraldton=(-28.778138, 114.615632)
)
loc_growthrates = map(locations) do loc
    map(f -> f[:H][Lat(Near(loc[1])), Lon(Near(loc[2]))], growth_output)
end
growthplot = plot(cat(loc_growthrates...); 
    labels=string.(keys(locations)), xguide="Time", yguide="Population"
)

init = (H=hostpop, P=parapop)
host_para_output = ArrayOutput(init; kwargs...)
sim!(host_para_output, localdisp, wind, host_para_growth)

p7 = plot(host_para_output[300][:H]; yguide="Host-Parasite", simplot_opts..., xtickfontcolor=RGB(1.0))
p8 = plot(host_para_output[600][:H]; yguide="Host-Parasite", simplot_opts..., xtickfontcolor=RGB(1.0))
hp_plots = plot(p7, p8; layout=(2, 1), size=(600, 600))

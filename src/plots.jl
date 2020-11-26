
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
    aux=(; rH=growthrates1, rP=growthrates1), 
    mask=parent(mask),
    tspan=tspan
)
wind_output = ArrayOutput(init; kwargs...)
localdisp_output = ArrayOutput(init; kwargs...)
combined_output = ArrayOutput(init; kwargs...)
sim!(wind_output, wind, growth)
sim!(localdisp_output, localdisp, growth)
sim!(combined_output, localdisp, wind, growth)

simplot_opts = (; 
    yguide_position=:right,
    ytickfontcolor=RGB(0.5), 
    clims=(0, carrycap),
    legend=:none, plot_opts...
)
p3 = plot(localdisp_output[600][:H]; yguide="Local", simplot_opts..., xtickfontcolor=RGB(1.0))
p4 = plot(wind_output[600][:H]; yguide="Wind", simplot_opts..., xtickfontcolor=RGB(1.0))
p5 = plot(combined_output[600][:H]; yguide="Combined", simplot_opts..., xtickfontcolor=RGB(0.5))
map(p -> plot!(p, shp; shape_opts...), (p3, p4, p5))
simplots = plot(p3, p4, p5; layout=(3, 1), size=(600, 900))
# savefig("output/simplots.png")

init = (H=hostpop, P=parapop)
host_para_output = ArrayOutput(init; kwargs...)
sim!(host_para_output, localdisp, wind, host_para_growth)

p6 = plot(host_para_output[300][:H]; yguide="Host-Parasite", simplot_opts..., xtickfontcolor=RGB(1.0))
p7 = plot(host_para_output[600][:H]; yguide="Host-Parasite", simplot_opts..., xtickfontcolor=RGB(1.0))
hp_plots = plot(p6, p7; layout=(2, 1), size=(600, 600))

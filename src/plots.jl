using Shapefile, Plots, Colors, ColorSchemes
using GLMakie
using Plots: px
Plots.theme(:default)

basedir = "/home/raf/julia/MethodsPaper"
include(joinpath(basedir, "src", "models.jl"))
include(joinpath(basedir, "src", "data.jl"))

# Plot setup
shp = Shapefile.Handle("/home/raf/.julia/dev/MethodsPaper/data/ausborder_polyline.shp")
shape_opts = (lw=0.5, color=RGB(0.3))
plot_opts = (; 
    showaxis=false,
    bottom_margin=0px, 
    foreground_color=RGB(0.3),
    tickfontsize=7, 
    title="", 
    xlims=(112, 155), ylims=(-44, -8)
)

#### Simulation plots ##########################################

tspan = DateTime(2020, 1):Day(7):DateTime(2030, 1)
kw = (; 
    aux=(; rH=rH, rP=rP), 
    mask=parent(mask),
    tspan=tspan
)
wind_output = ArrayOutput(initdata; kw...)
localdisp_output = ArrayOutput(initdata; kw...)
combined_output = ArrayOutput(initdata; kw...)
growth_output = ArrayOutput((; H=copy(hostinit) .= carrycap .* mask); kw...)
host_para_output = ArrayOutput(initdata; kw...)
sim!(wind_output, wind, allee, growth)
sim!(localdisp_output, localdisp, allee, growth)
sim!(combined_output, localdisp, wind, allee, growth)
sim!(growth_output, allee, growth)
sim!(host_para_output, localdisp, wind, allee, growth, localdispP, alleeP, parasitism)

simplot_opts = (; 
    # left_margin=-50px,
    clims=(0, carrycap),
    legend=:none,
    showaxis=false,
    ticks=false,
    # xtickfontcolor=RGB(0.5),
    # ytickfontcolor=RGB(0.5),
    yguide="",
    plot_opts...
)
t = Ti(Near(DateTime(2022, 2)))
p3 = Plots.plot(localdisp_output[t][:H]; xguide="Local", simplot_opts...)
p4 = Plots.plot(wind_output[t][:H]; xguide="Wind", simplot_opts...)
p5 = Plots.plot(combined_output[t][:H]; xguide="Combined", simplot_opts...)
map(p -> Plots.plot!(p, shp; shape_opts...), (p3, p4, p5))
simplots = Plots.plot(p3, p4, p5; layout=(1, 3), size=(900, 300))

savefig("output/simplots.png")

#### Host Parasite simulation plots #############################

pyplot()
t1 = DateTime(2021, 1)
t2 = DateTime(2025, 1)
t3 = DateTime(2030, 1)
hp_plot_opts = (; 
    ytickfontcolor=RGB(1.0),
    simplot_opts...
)
col_p = cgrad(ColorSchemes.diverging_linear_bjr_30_55_c53_n256)
no_x = (xguide="", xtickfontcolor=RGB(1.0))#, bottom_margin=-30px) 
hp1, hp2, hp3 = map((t1, t2, t3), (no_x, no_x, ())) do t, kw
    Plots.plot(host_para_output[Ti(Near(t))][:H]; xguide="Host", hp_plot_opts..., kw..., yguide=year(t)),
    Plots.plot(host_para_output[Ti(Near(t))][:P]; color=col_p, xguide="Parasite", hp_plot_opts..., ytickfontcolor=RGB(1.0), kw...)
    #plot(combined_output[Ti(Near(t))][:H]; xguide="Host only", hp_plot_opts..., ytickfontcolor=RGB(1.0), kw...)
end
map(p -> Plots.plot!(p, shp; shape_opts...), (hp1..., hp2..., hp3...))
host_para_plot = Plots.plot(hp1..., hp2..., hp3...; layout=(3, 2), size=(900, 900))

savefig("output/hostparaplots.png")


#### Growth line plot for locations ############################

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
plotspan = Between(DateTime(2020), DateTime(2021))
loc_populations = map(locations) do loc
    coords = Lat(Contains(loc[1])), Lon(Contains(loc[2]))
    map(f -> f[:H][coords...] + 1, growth_output[plotspan])
end
# plot(growth_output[end][:H])
frames = map(f -> f[:H], growth_output[Ti(plotspan)])
# grid = get(scheme, cat(frames...; dims=Ti(1:length(frames))) ./ carrycap)
grid = cat(frames...; dims=dims(frames, Ti)) ./ carrycap
scene = GLMakie.volume(parent(grid); algorithm=:absorptionrgba)
save("iso.png", scene)
Axis(scene)

function dimpoints(g, dim) 
    dnum = dimnum(g, dim) 
    ind = index(g, dnum)
    [ind[Tuple(I)[dnum]] for I in CartesianIndices(g) if g[I] != 0] 
end

x = dimpoints(grid, Lon())
y = dimpoints(grid, Lat())
z = dimpoints(grid, Ti())
Plots.plot(grid[Ti(10)])
colors = [c for c in grid if c != 0]
GLMakie.meshscatter(x, y, z)


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

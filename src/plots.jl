using Shapefile, Plots, Colors, ColorSchemes, Pkg
using Plots: px
Plots.theme(:default)

# The directory for this project
basedir = dirname(Pkg.project().path)

# Inlude models and data scripts
include(joinpath(basedir, "src", "rules.jl"))
include(joinpath(basedir, "src", "data.jl"))

# Define rules and simulation data
growth, allee, localdisp, wind, parasitism, localdisp_p, allee_p, carrycap = define_rules()
init_h, init_p, rate_h, rate_p, mask = load_data(basedir, carrycap)

# Define NamedTuples of initialisation and auxilary data for the simulations.
initdata = (; H=init_h, P=init_p)
auxdata = (; rH=rate_h, rP=rate_p)

# Save gifs, as in the paper but with colors
output = GifOutput(initdata; 
    filename=joinpath(basedir, "output", "host_parasitoid.gif"), 
    tspan=DateTime(2020, 1):Day(7):DateTime(2030, 1),
    aux=(rH=rate_h, rP=rate_p,), mask=mask, # See S3
    minval=(0.0, 0.0), maxval=(carrycap, carrycap),
    scheme=(ColorSchemes.inferno, ColorSchemes.diverging_linear_bjr_30_55_c53_n256),
    text=nothing,
) 
sim!(output, (wind, localdisp, allee, growth, localdisp_p, allee_p, parasitism))

# Plot setup
shp = Shapefile.Handle(joinpath(basedir, "data", "ausborder_polyline.shp"))
shape_opts = (lw=0.5, color=RGB(0.3))
plot_opts = (; 
    showaxis=false,
    bottom_margin=0px, 
    foreground_color=RGB(0.3),
    tickfontsize=7, 
    title="", 
    xlims=(112, 155), ylims=(-44, -8)
)

#### Host simulation plots ##########################################
tspan = DateTime(2020, 1):Day(7):DateTime(2030, 1)
kw = (; 
    aux=auxdata, mask=mask, tspan=tspan
)
# Define outputs
wind_output = ArrayOutput(initdata; kw...)
localdisp_output = ArrayOutput(initdata; kw...)
combined_output = ArrayOutput(initdata; kw...)
host_para_output = ArrayOutput(initdata; kw...)
# Run simulations
sim!(wind_output, wind, allee, growth)
sim!(localdisp_output, localdisp, allee, growth)
sim!(combined_output, localdisp, wind, allee, growth)
sim!(host_para_output, localdisp, wind, allee, growth, localdisp_p, allee_p, parasitism)

# Default map options
simplot_opts = (; 
    color=:inferno,
    clims=(0, carrycap),
    legend=:none,
    showaxis=false,
    ticks=false,
    yguide="",
    plot_opts...
)

# Make simulation plots
t = Ti(Near(DateTime(2022, 2)))
p3 = Plots.plot(localdisp_output[t][:H]; xguide="Local", simplot_opts...)
p4 = Plots.plot(wind_output[t][:H]; xguide="Wind", simplot_opts...)
p5 = Plots.plot(combined_output[t][:H]; xguide="Combined", simplot_opts...)
map(p -> Plots.plot!(p, shp; shape_opts...), (p3, p4, p5))
simplots = Plots.plot(p3, p4, p5; layout=(1, 3), size=(900, 300))

savefig(joinpath(basedir, "output", "simplots.png"))

#### Host/Parasitoid simulation plots #############################

t1 = DateTime(2021, 1)
t2 = DateTime(2025, 1)
t3 = DateTime(2030, 1)
hp_plot_opts = (; 
    ytickfontcolor=RGB(1.0),
    simplot_opts...
)
col_p = cgrad(ColorSchemes.diverging_linear_bjr_30_55_c53_n256)
no_x = (xguide="", xtickfontcolor=RGB(1.0))
hp1, hp2, hp3 = map((t1, t2, t3), (no_x, no_x, ())) do t, kw
    hp = Plots.plot(host_para_output[Ti(Near(t))][:H]; xguide="Host", hp_plot_opts..., kw..., yguide=year(t)),
         Plots.plot(host_para_output[Ti(Near(t))][:P]; xguide="Parasitoid", hp_plot_opts..., ytickfontcolor=RGB(1.0), kw..., color=:cividis)
    map(x -> Plots.plot!(x, shp; shape_opts...), hp)
end
host_para_plot = Plots.plot(hp1..., hp2..., hp3...; layout=(3, 2), size=(900, 900))

savefig(joinpath(basedir, "output", "hostparaplots.png"))

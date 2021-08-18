using ColorSchemes, Pkg, DynamicGridsGtk

# Warning: Due to a threading bug in GTK, do not run benchmarks or use
# proc=ThreadedCPU after running this file. They will crawl!

# The directory for this project
basedir = dirname(Pkg.project().path)

# Inlude models and data scripts
include(joinpath(basedir, "src", "rules.jl"))
include(joinpath(basedir, "src", "data.jl"))

# Define rules and simulation data
growth, allee, localdisp, wind, parasitism, localdisp_p, allee_p, carrycap = define_rules()
init_h, init_p, rate_h, rate_p, mask = load_data(basedir, carrycap)

# Define NamedTuples auxilary data for the simulations.
auxdata = (; rH=rate_h, rP=rate_p)

tspan = DateTime(2020, 1):Day(7):DateTime(2040, 1)

viewer = CopyTo{:rH}(Aux(:rH))
host_ruleset = Ruleset(growth, localdisp, allee, viewer)

output = GtkOutput((H=init_h, rH=zero(init_h)); 
    ruleset=host_ruleset,
    aux=auxdata, mask=mask,
    tspan=tspan, 
    store=false,
    fps=100,
    minval=(0.0f0, 0.0f0),
    maxval=(carrycap, maximum(rate_h)),
    scheme=(ColorSchemes.autumn, ColorSchemes.inferno),
)
sim!(output, (growth, wind, localdisp, allee, viewer)) 

host_para_ruleset = Ruleset(wind, localdisp, growth, allee, localdisp_p, allee_p, parasitism)

output = GtkOutput((; H=init_h, P=init_p),
    filename="output/hostparasite.gif", 
    ruleset=host_para_ruleset,
    aux=auxdata, mask=mask,
    tspan=tspan,
    fps=200,
    store=true,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, carrycap),
    scheme=(ColorSchemes.autumn, ColorSchemes.jet),
)
sim!(output, host_para_ruleset)

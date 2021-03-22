using ColorSchemes, DynamicGridsGtk#, DynamicGridsInteract

# Warning: Due to a threading bug in GTK, do not run benchmarks or use
# proc=ThreadedCPU after running this file. They will crawl!

# The directory for this project
basedir = realpath(@__DIR__)
# Inlude models and data scripts
include(joinpath(basedir, "src", "models.jl"))
include(joinpath(basedir, "src", "data.jl"))

tspan = DateTime(2020, 1):Day(7):DateTime(2040, 1)

viewer = CopyTo{:rH}(Aux(:rH))
host_ruleset = Ruleset(growth, localdisp, allee, viewer)

# Comment/uncomment these to choose the interface
# To generate sliders in ElectronOuput, model parameters
# need to be defined as: `Param(val; bounds=(lower, upper))` - see the 
# ModelParameters.jl docs for details.

# output = ElectronOutput((H=init_h, rH=zero(init_h)); 
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

# output = ElectronOutput(initdata;
output = GtkOutput((H=init_h, P=init_p);
    filename="output/hostparasite.gif", 
    ruleset=host_para_ruleset,
    aux=auxdata, mask=mask,
    tspan=tspan,
    fps=100,
    store=true,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, carrycap),
    scheme=(ColorSchemes.autumn, ColorSchemes.jet),
)
sim!(output, ruleset)

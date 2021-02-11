
# Warning: Due to a threading bug in GTK, do not run benchmarks or use
# ThreadedCPU during or after running this file. They will crawl.
using ColorSchemes, DynamicGridsGtk

tspan = DateTime(2020, 1):Day(7):DateTime(2040, 1)

viewer = CopyTo{:rH}(Aux(:rH))
ruleset = Ruleset(growth, localdisp, allee, viewer)
output = GtkOutput((H=hostinit, rH=zero(hostinit)); 
    ruleset=ruleset,
    aux=auxdata,
    mask=parent(mask),
    tspan=tspan, 
    store=false,
    fps=100,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, maximum(rH)),
    scheme=(ColorSchemes.autumn, ColorSchemes.inferno),
)
sim!(output, (growth, wind, localdisp, allee, viewer)) 

ruleset = Ruleset(wind, localdisp, growth, allee, localdispP, alleeP, parasitism)
# ruleset = Ruleset(randomgrid, gpu_wind, localdisp, growth, allee, localdispP, alleeP, parasitism; proc=CuGPU())

output = GtkOutput(initdata;
# output = ElectronOutput((H=hostpop, P=parapop);#, rH=zero(rH[Ti(1)])); 
    filename="output/hostparasite.gif", 
    ruleset=ruleset,
    aux=auxdata,
    mask=parent(mask),
    tspan=tspan,
    fps=1000,
    store=true,
    minval=(0.0f0, 0.0f0, nothing), 
    maxval=(carrycap, carrycap, nothing),
    scheme=(ColorSchemes.autumn, ColorSchemes.jet, Greyscale()),
)
sim!(output, ruleset)
# savegif("output/hostparasite.gif", output)

using DynamicGridsInteract, ColorSchemes
using DynamicGridsGtk

viewer = CopyTo{:rH}(Aux(:rH))
ruleset = Ruleset(growth, allee, viewer)
# output = GtkOutput((H=hostpop, rH=zero(hostpop)); 
output = ElectronOutput((H=fill(1f3, size(hostpop)), rH=zero(hostpop)); 
    ruleset=ruleset,
    aux=(; rH=rH),
    mask=parent(mask),
    tspan=tspan, 
    store=true,
    fps=50,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, maximum(rH)),
    scheme=(ColorSchemes.autumn, ColorSchemes.inferno),
)

sim!(output, (growth, viewer)) 

ruleset = Ruleset(wind, localdisp, localdispP, allee, alleeP, host_parasite)
# output = GtkOutput((H=hostpop, P=parapop, rH=zero(rH[Ti(1)])); 
output = ElectronOutput((H=hostpop, P=parapop);#, rH=zero(rH[Ti(1)])); 
    ruleset=ruleset,
    aux=(; rH=rH, rP=rH),
    mask=parent(mask),
    tspan=tspan,
    fps=100,
    store=true,
    minval=(0.0f0, 0.0f0), 
    maxval=(carrycap, carrycap),
    scheme=(ColorSchemes.autumn, ColorSchemes.autumn),
)

sim!(output, ruleset)

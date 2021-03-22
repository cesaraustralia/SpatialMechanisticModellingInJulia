Scripts to accompany the paper _"A working guide to spatial mechanistic modelling
in Julia"_ by Rafael Schouten, Virgile Baudrot, Paul Umina and James Maino.

To build all the plots from the command line, run from the repo directory:

```bash
julia --project=. src/plots.jl
```

Otherwise they can be stepped through in Atom or VSCode, which may be better
a better way to run the simulations defined in outputs.jl. In this case be sure 
to do:

```julia
] 
activate .
instantiate
```

To activate the project before running anything.

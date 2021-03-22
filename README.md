Scripts to accompany the paper _"A working guide to spatial mechanistic modelling
in Julia"_ by Rafael Schouten, Virgile Baudrot, Paul Umina and James Maino.

To build the growth-rate maps, run from the directory this file is in:

```bash
julia --project=. src/growthmaps.jl
```

To build all the plots:

```bash
julia --project=. src/plots.jl
```

Otherwise scripts can be stepped through in Atom or VSCode, which may be better
a better way to run the simulations defined in outputs.jl. In this case be sure 
to do:

```julia
] 
activate .
instantiate
```

To activate the project before running anything.

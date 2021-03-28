Scripts to accompany the paper _"A working guide to spatial mechanistic
modelling in Julia"_ by Rafael Schouten, James Maino, Virgile Baudrot and Paul Umina.


## Setup

[Download](https://julialang.org/downloads/) and add Julia to your system.

We recommend Julia 1.6 or higher. 

Then, download or clone this repository from github.

## Use

To e.g. build all the plots from the command line, navigate to the repository folder
and run:

```bash
julia --project=. src/plots.jl
```

Otherwise the scripts can be stepped through in Atom with Juno or VSCode with
the julia-vscode plugin. This may be a better way to run the simulations
defined in outputs.jl, and explore parameterisations. 

Be sure to start Julia in the project folder. Get a REPL using `Ctrl+shift+p` and search
for julia repl in vscode, 

and do:

```julia
] 
activate .
instantiate
```

Or use:

```julia
julia --project=.
```

To activate the project before running any of the code. This will ensure that
you have the exact versions of dependencies, which means these scripts will
continue to work long into the future.

## Included scripts:

- models.jl: defines the rule components for the models we simulate.
- data.jl: prepares simulation data.
- plots.jl: creates most plots seen in paper from models.jl/data.jl.
- benchmarks.jl: runs simulation bechmarks and creates benchmark plot.
- outputs.jl: provides some interface visualisations of the simulations.
- growthmaps.jl: shows how growth-rate data used to drive growth models is generated.


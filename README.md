Scripts to accompany the paper _"A working guide to spatial mechanistic
modelling in Julia"_ by Rafael Schouten, James Maino, Virgile Baudrot and Paul
Umina.

These scripts should be reproducible on all systems, including all the binaries
used by various packages.

## Setup

[Download](https://julialang.org/downloads/) and add Julia to your system.

Julia 1, from v1.6 or higher is required.

Then, download or clone this code repository locally. To install exact versions
of all packages, start a julia REPL in the downloaded project folder, and in
package mode (with `]`) activate and instantiate the project:

```julia
] 
activate .
instantiate
```

Or from a command line, in the repository directory:

```bash
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
```

## Use

To e.g. build all the plots from the command line, in the repository directory
you can run:

```bash
julia --project=. src/plots.jl
```

Scripts can also be stepped through in an editor like VSCode, with the
julia-vscode plugin. This may be a better way to explore and analyse the code,
and to run the simulations in `outputs.jl`. 

Be sure to start Julia in the project folder. The project directory should be
activated automatically in VSCode, if not, or in other editors, in the REPL
type:

```julia
] 
activate .
```

This will ensure that you have the exact versions of all dependencies, which
means these scripts will continue to work into the future.

## Included scripts

To run directly:
- plots.jl: creates most plots seen in paper from rules.jl/data.jl.
- benchmarks.jl: runs simulation bechmarks and creates benchmark plot.
- outputs.jl: provides some interface visualisations of the simulations.
- growthmaps.jl: builds the growth-rate data used to drive growth models.

Called by the other scripts:
- rules.jl: defines the rules used in the simulations.
- data.jl: prepares simulation data.

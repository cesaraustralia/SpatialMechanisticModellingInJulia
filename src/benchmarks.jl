using BenchmarkTools, Plots, CUDAKernels, CUDA, Random, Pkg, Adapt
using DynamicGrids: CuGPU, CPU, GPU, SimData
using Plots: px
theme(:vibrant)
CUDA.allowscalar(false)
# The directory for this project
basedir = dirname(Pkg.project().path)
# Inlude models and data scripts
include(joinpath(basedir, "src", "rules.jl"))
include(joinpath(basedir, "src", "data.jl"))

##### GPU wind rule #####

# We can't use `rand` directly _inside_ a GPU kernel, as of February 2021.
# To get around this shortcoming we add a grid to hold random numbers and
# call the random generator from outside a kernel using a `SetGrid` rule that
# calls a function on the entire grid. All the work is still done on the GPU.
gpu_wind = SetNeighbors{Tuple{:rand,:H}}() do data, hood, (rand, H), I
    # Ignore empty cells
    H > zero(H) || return nothing
    # Randomise a dispersal event using cell value of rand
    rand < 0.02f0 || return nothing
    # Get random values for jump from the surrounding grid
    randi, randj = (I[1], I[2] + 1), (I[1] + 1, I[2] + 1)
    dest = I .+ (unsafe_trunc(Int, data[:rand][randi...] * 40 - 20),
                 unsafe_trunc(Int, data[:rand][randj...] * 40 - 20))
    # Update spotted cell if it's on the grid
    if isinbounds(data[:H], dest)
        @inbounds add!(data[:H], H / 10, dest...)
        @inbounds sub!(data[:H], H / 10, I...)
    end
    nothing
end

_rand!(A) = rand!(A)
_rand!(A::CuArray) = CUDA.rand!(A)

# A rule that just fills the random grid
randomgrid = SetGrid{Tuple{},:rand}() do randgrid
    # `parent` is so we write rand to the whole grid, including any padding
    _rand!(parent(randgrid))
end


##### Benchmark setup methods #####

function setupsim(rules;
    tspan=DateTime(2020, 1):Week(1):DateTime(2022, 1),
    opt=NoOpt(), proc=SingleCPU(), size_ag,
    growthmin=-5.0f0, growthmax=0.2f0,
    output_type=ResultOutput,
)
    # Take top right corner to make a square - but most of australia
    ax = 1:200, lastindex(init_h, 2)-199:lastindex(init_h, 2)
    sze = size_ag[1]
    ag = size_ag[2]
    hpop = ag(init_h[ax...])
    @assert size(hpop, 1) == size(hpop, 2) == sze
    init = (;
        H=ag(init_h[ax...]),
        P=ag(init_p[ax...]),
        rand=rand(Float32, sze, sze)
    )
    o = output_type(init;
        aux=(; rH=ag(rate_h[ax..., :]), rP=ag(rate_p[ax..., :])),
        mask=Array(ag(mask[ax...])),
        tspan=tspan,
        style=Braile(), color=:blue, fps=100,
    )
    if proc isa GPU
        opt isa SparseOpt && error("Sparse GPU")
        proc_rules = rules[:gpu]
        # Copy data to GPU so it's not part of the timing - otherwise
        # we are just benchmarking CPU/GPU data transfer.
        o = adapt(CuArray, o)
        rs = Ruleset(proc_rules; opt=opt, proc=proc)
        sd = adapt(CuArray, SimData(o, rs))
    else
        proc_rules = rules[:cpu]
        rs = Ruleset(proc_rules; opt=opt, proc=proc)
        sd = SimData(o, rs)
    end
    return (; o, rs, sd)
end
function suite(rules)
    suite = BenchmarkGroup()
    for opt in keys(opts)
        suite[opt] = BenchmarkGroup([:proc, :size])
        for proc in keys(procs), s in sizes
            # SparseOpt is not implemented for GPU yet
            proc == :gpu && opt == :sparseopt && continue
            suite[opt][proc, first(s)] = @benchmarkable sim!(o, rs; simdata=sd) setup=((o, rs, sd)=setupsim($rules; size_ag=$s, proc=$(procs[proc]), opt=$(opts[opt])))
        end
    end
    suite
end


##### Scenarios and benchmark #####

procs = (single=SingleCPU(), threaded=ThreadedCPU(), gpu=CuGPU(),)
opts = (noopt=NoOpt(), sparseopt=SparseOpt())
# We use a an anonymous function to scale a 200*200 array to the specified size
sizes = (
    100 => A -> aggregate(Center(), A, (Lat(2), Lon(2), Ti(1))),
    200 => A -> A,
    400 => A -> disaggregate(Center(), A, (Lat(2), Lon(2), Ti(1))),
    800 => A -> disaggregate(Center(), A, (Lat(4), Lon(4), Ti(1))),
    1600 => A -> disaggregate(Center(), A, (Lat(8), Lon(8), Ti(1))),
)

# CPU/GPU rules are separated as the best optimisations are not the same on CPU and GPU.
# Chain has better effects on GPU than on GPU, and GPU needs random numbers generated
# in a GridRule on a separate grid, due to the current lack of a `rand` function
# inside GPU kernels.
rulegroups = (
    Wind=(cpu=(wind, allee, growth), gpu=(randomgrid, gpu_wind, allee, growth)),
    Local=(cpu=(localdisp, allee, growth), gpu=(localdisp, allee, growth)),
    Combined=(cpu=(wind, localdisp, allee, growth), gpu=(randomgrid, gpu_wind, localdisp, allee, growth)),
    Parasitism=(cpu=(wind, localdisp, growth, allee, localdisp_p, allee_p, parasitism),
                gpu=(randomgrid, gpu_wind, localdisp, localdisp_p, Chain(allee, growth, allee_p, parasitism))),
)

suites = map(suite, rulegroups)
map(tune!, suites)
results = map(su -> run(su, verbose=true), suites)

##### Benchmark plot #####

function plotbench(b, key)
    @show b
    p = plot()
    i = 2
    for prockey in Tuple(keys(procs))
        for optkey in Tuple(keys(opts))
            if optkey == :sparseopt
                prockey == :gpu && continue
                lab = uppercasefirst(string("Sparse ", prockey))
            else
                lab = if prockey == :gpu
                    uppercase(string(prockey))
                else
                    uppercasefirst(string(prockey))
                end
            end
            times = map(sizes) do (s, ag)
                median(b.data[string(optkey)].data[(string(prockey), s)].times) / 1e9
            end
            plot!(p, collect(map(first, sizes)), collect(times);
                color=i, label=lab, pallete=:Dark2_5, opacity=0.6
            )
            i += 1
        end
    end
    plot(p;
        title=string(key),
        ylims=(1e-3, 10.0),
        yaxis=:log,
        titlefontsize=10,
        guidefontsize=9,
        tickfontsize=8,
        linewidth=2,
        legend=(key == :Parasitism ? :bottomright : :none),
        tickfontcolor=RGB(0.2),
        legendfontcolor=RGB(0.2),
        minorgrid=false,
        grid=true,
        xrotation=60,
        xticks=[100, 200, 400, 800, 1600],
        xlabel="Grid size",
        ylabel=(key in (:Wind, :Combined) ? "Time in seconds" : ""),
    )
end
plot(map(plotbench, results, keys(results))...; layout=(2, 2), size=(530, 600))

savefig("output/benchmarks.png")

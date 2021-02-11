using BenchmarkTools, Plots, KernelAbstractions, CUDA, Random, Adapt
using DynamicGrids: CuGPU, CPU, GPU

# CUDA tweaks
CUDA.allowscalar(false)
# Use CUDA rand
_rand!(A::CuArray) = CUDA.rand!(A)

basedir = "/home/raf/julia/MethodsPaper"
include(joinpath(basedir, "src", "models.jl"))
include(joinpath(basedir, "src", "data.jl"))


function setupsim(rules;
    tspan=DateTime(2020, 1, 1):Day(1):DateTime(2020, 4, 10),
    opt=NoOpt(), proc=SingleCPU(), size_ag,
    growthmin=-5.0f0, growthmax=0.2f0,
    output_type=ResultOutput
)
    # Take top right corner to make a square - but most of australia
    ax = 1:200, lastindex(hostpop, 2)-199:lastindex(hostpop, 2)
    sze = size_ag[2]
    ag = size_ag[2]
    hpop = parent(ag(hostpop[ax...]))
    @assert size(hpop, 1) == size(hpop, 2) == sze
    init = (; 
        H=parent(ag(hostpop[ax...])), 
        P=parent(ag(parapop[ax...])),
    )
    o = output_type(init, rand=rand(Float32, sze, sze); 
        aux=(; rH=ag(rH[ax..., :]), rP=ag(rP[ax..., :])), 
        mask=Array(ag(mask[ax...])),
        tspan=tspan,
        style=Braile(), color=:blue, fps=100,
    )
    proc_rules = rules[proc isa CPU ? :cpu : :gpu]
    rs = Ruleset(proc_rules; opt=opt, proc=proc)
    return (; o, rs)
end

function suite(rules)
    suite = BenchmarkGroup()
    for opt in keys(opts)
        suite[opt] = BenchmarkGroup(["proc", "size"])
        for proc in keys(procs), s in sizes
            # SparseOpt is not implemented for GPU yet
            proc isa CuGPU && opt isa SparseOpt && continue
            suite[opt][proc, s] = @benchmarkable sim!(o, rs) setup=((o, rs)=setupsim($rules; size=$s, proc=$(procs[proc]), opt=$(opts[opt])))
        end
    end
    suite
end

procs = (single=SingleCPU(), threaded=ThreadedCPU(), gpu=CuGPU())
opts = (noopt=NoOpt(), sparseopt=SparseOpt())
sizes = [
    100 => A -> aggregate(Center(), A, Lat(2), Lon(2)), 
    200 => A -> A,  
    400 => A -> disaggregate(Center(), A, Lat(2), Lon(2)), 
    800 => A -> disaggregate(Center(), A, Lat(4), Lon(4)), 
]

# CPU/GPU rules 
# These are separate as the best optimisations are not the same on CPU and GPU. 
# Chain has better effects on GPU than on GPU, and GPU needs random numbers generated 
# in a GridRule # on a separate grid, due to the current lack of a `rand` function 
# inside GPU kernels.
rulegroups = (
    Wind=(cpu=(wind, allee, growth), gpu=(randomgrid, gpu_wind, allee, growth)), 
    Local=(cpu=(localdisp, allee, growth), gpu=(localdisp, allee, growth)), 
    Combined=(cpu=(wind, localdisp, allee, growth), gpu=(randomgrid, gpu_wind, localdisp, allee, growth)),
    Parasitism=(cpu=(wind, localdisp, growth, allee, localdispP, alleeP, parasitism), 
                gpu=(randomgrid, gpu_wind, localdisp, localdispP, Chain(allee, growth, alleeP, parasitism))),
)

suites = map(suite, rulegroups)
map(tune!, suites)
results = map(s -> run(s, verbose=true), suites)

# Run all variants in the REPL to see that they make sense
for key in keys(rulegroups), size in sizes, opt in values(opts), proc in values(procs)
    @show key
    o, rs = setupsim(rulegroups[key], size=size)#, output_type=REPLOutput)
    sim!(o, rs; opt=opt, proc=proc)
end

@time results = map(bench, rulegroups)
results[:Combined][@tagged 400]
# minimum(results[:Combined].data["sparseopt"].data[("gpu", 400)].times) / 1e9

# tspan = DateTime(2020, 1):Week(1):DateTime(2020, 1) + Week(100)
# output = ArrayOutput((; H=parent(hostpop), P=parent(parapop), rand=rand(Float32, size(hostpop))); 
#     aux=(; rH=rH, rP=rP), 
#     mask=Array(mask),
#     tspan=tspan
# )

# sd = DynamicGrids.SimData(output, Ruleset(rulegroups[:Wind]))

# Flatten.fieldnameflatten(sd, Union{Array,BitArray}, SArray)
# Flatten.flatten(sd, Union{Array,BitArray}, Union{SArray,Dict})
# using Adapt, Flatten, StaticArrays

function plotbench(b, key)
    p = plot()
    i = 2
    for prockey in Tuple(keys(procs))
        for optkey in Tuple(keys(opts))
            if optkey == :sparseopt
                prockey == :gpu && continue
                lab = uppercasefirst(string("Sparse ", prockey))
            else
                lab = uppercasefirst(string(prockey))
            end
            times = map(sizes) do s
                median(b.data[optkey].data[(prockey, s)].times) / 1e9
            end
            plot!(p, collect(sizes), collect(times); 
                color=i, label=lab, pallete=:Dark2_5, opacity=0.6
            )
            i += 1
        end
    end
    plot(p; 
        title=string(key), 
        ylims=(1e-3, 1.0), 
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
        xlabel=(key in (:Combined, :Parasitism) ? "Size" : ""),
        # xformatter=(key in (:Combined, :Parasitism) ? _->"" : false),
        ylabel=(key in (:Wind, :Combined) ? "Time in seconds" : ""),
        # yticks=(key in (:Wind, :Combined) ? true : false),
    )
end
# pyplot()
theme(:vibrant)
using Plots: px
plot(map(plotbench, results, keys(results))...; layout=(2, 2), size=(500, 500))

savefig("output/benchmarks.png")


### Other ########################################

# using ConstructionBase
# ConstructionBase.constructorof(::Type{<:LinRange}) = _linrange
# _linrange(start, stop, len, lendiv) = LinRange(start, stop, len)

# hostpop_disag, parapop_disag, rH_disag, rP_disag, mask_disag = 
#     map((hostpop, parapop, rH, rP, mask)) do A
#         A = set(A, Lat=NoIndex(), Lon=NoIndex())
#         # A = hasdim(A, Ti) ? set(A, Ti=NoIndex()) : A
#         GeoData.disaggregate(Center(), A, (Lat(4), Lon(4), Ti(1)))
#     end
# A = set(rH, Lat=NoIndex(), Lon=NoIndex(), Ti=NoIndex())
# output = ArrayOutput((; H=hostpop_disag, P=parapop_disag, rand=rand(Float32, size(hostpop_disag))); 
#     aux=(; rH=rH_disag, rP=rP_disag), 
#     mask=Array(mask_disag), 
#     tspan=tspan
# )

tspan = DateTime(2020, 1):Week(1):DateTime(2020, 1) + Week(100)
output = ArrayOutput(initdata; aux=auxdata, mask=Array(mask), tspan=tspan)

@time sim!(output, wind);
@time sim!(output, wind; opt=SparseOpt());

@time sim!(output, growth);
@time sim!(output, growth; proc=ThreadedCPU());
@time sim!(output, growth; proc=CuGPU());
@time sim!(output, localdisp);
@time sim!(output, localdisp; proc=ThreadedCPU());
@time sim!(output, localdisp; proc=CuGPU());
@time sim!(output, (localdisp, growth, allee));
@time sim!(output, (localdisp, growth, allee); proc=ThreadedCPU());
@time sim!(output, Chain(localdisp, growth, allee));
@time sim!(output, Chain(localdisp, growth, allee); proc=ThreadedCPU());
@time sim!(output, (localdisp, growth, allee); proc=CuGPU());
@time sim!(output, Chain(localdisp, growth, allee); proc=CuGPU());
@time sim!(output, randomgrid, gpu_wind, growth);

@time sim!(output, wind, Chain(growth, allee));
@time sim!(output, randomgrid, gpu_wind, growth; proc=ThreadedCPU());
@time sim!(output, randomgrid, gpu_wind, growth; proc=CuGPU());
@time sim!(output, randomgrid, gpu_wind, Chain(growth, allee); proc=ThreadedCPU());
@time sim!(output, randomgrid, gpu_wind, Chain(growth, allee); proc=CuGPU());

@time sim!(output, randomgrid, gpu_wind, Chain(localdisp, growth, allee));
@time sim!(output, randomgrid, gpu_wind, Chain(localdisp, growth, allee); proc=ThreadedCPU());
@time sim!(output, randomgrid, gpu_wind, Chain(localdisp, growth, allee); proc=CuGPU());
@time sim!(output, randomgrid, gpu_wind, localdisp, growth, allee);
@time sim!(output, randomgrid, gpu_wind, localdisp, growth, allee; proc=ThreadedCPU());
@time sim!(output, randomgrid, gpu_wind, localdisp, growth, allee; proc=CuGPU());

@time sim!(output, growth, opt=SparseOpt());
@time sim!(output, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, growth; proc=CuGPU(), opt=SparseOpt());
@time sim!(output, localdisp, growth, allee, opt=SparseOpt());
@time sim!(output, Chain(localdisp, growth, allee), opt=NoOpt());
@time sim!(output, Chain(localdisp, growth, allee), opt=SparseOpt());
@time sim!(output, localdisp, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, localdisp, growth; proc=CuGPU(), opt=SparseOpt());
@time sim!(output, randomgrid, gpu_wind, growth, opt=SparseOpt());
@time sim!(output, randomgrid, gpu_wind, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, randomgrid, gpu_wind, growth; proc=CuGPU(), opt=SparseOpt());
@time sim!(output, localdisp, randomgrid, gpu_wind, growth, opt=SparseOpt());
@time sim!(output, localdisp, randomgrid, gpu_wind, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, localdisp, randomgrid, gpu_wind, growth; proc=CuGPU(), opt=SparseOpt());

@time sim!(output, randomgrid, gpu_wind, localdisp, allee, localdispP, alleeP, parasitism; proc=SingleCPU());
@time sim!(output, randomgrid, gpu_wind, localdisp, allee, localdispP, alleeP, parasitism; proc=ThreadedCPU());
@time sim!(output, randomgrid, gpu_wind, localdisp, allee, localdispP, alleeP, parasitism; proc=CuGPU());
@time sim!(output, randomgrid, gpu_wind, Chain(localdisp, allee), Chain(localdispP, alleeP), parasitism; proc=SingleCPU());
@time sim!(output, randomgrid, gpu_wind, Chain(localdisp, allee), Chain(localdispP, alleeP), parasitism; proc=ThreadedCPU());
# @time sim!(output, randomgrid, gpu_wind, Chain(localdisp, allee), Chain(localdispP, alleeP), parasitism; proc=CuGPU());

ruleset = (randomgrid, gpu_wind, localdisp, para_localdisp, host_para_growth)
@time sim!(output, ruleset);
@time sim!(output, ruleset; proc=ThreadedCPU());
@time sim!(output, ruleset; proc=CuGPU());

using BenchmarkTools, Cthulhu

rule = Chain(localdisp, growth, allee)
rule = growth
tspan = DateTime(2020, 1):Week(1):DateTime(2020, 1) + Week(100)
output = ArrayOutput((; H=parent(hostpop), P=parent(parapop), rand=rand(Float32, size(hostpop))); 
    aux=(; rH=rH, rP=rP), mask=Array(mask), tspan=tspan
)
sd = DynamicGrids._initdata!(nothing, output.extent, Ruleset(rule))
sd = DynamicGrids._updatetime(sd, 1) |> s -> DynamicGrids.precalcrules(s, rules(s))
rule = first(rules(sd))
grid = sd[DynamicGrids.neighborhoodkey(rule)]
r = max(1, DynamicGrids.radius(rule))
T = eltype(grid)
S = 2r + 1
buffer = SArray{Tuple{S,S},T,2,S^2}(Tuple(zero(T) for i in 1:S^2))
@btime DynamicGrids._setbuffer($rule, $buffer)
DynamicGrids.maprule!(sd, rule)
@code_warntype DynamicGrids.applyrule(sd, rule, (H=1.0f0, P=1.0f0), (1,2))
@code_native DynamicGrids.applyrule(sd, rule, (H=1.0f0, P=1.0f0), (1,2))
@descend DynamicGrids.applyrule(sd, rule, (H=1.0f0, P=1.0f0), (1, 2))

function perf(f, args)
    pid = getpid()
    cmd = `perf $args --pid=$pid`
    proc = run(pipeline(cmd, stdout=stdout, stderr=stderr); wait=false)
    try
        return f()
    finally
        flush(stdout)
        flush(stderr)
        kill(proc, Base.SIGINT)
        wait(proc)
    end
end

perf(f, `stat -e L1-dcache-load-misses,L1-dcache-loads,L1-dcache-stores,L1-icache-load-misses`)
perf(f, `stat -e LLC-load-misses,LLC-loads,LLC-stores`)

perf(f, ``) 

f() = (@btime sim!($output, $localdisp, $growth; proc=ThreadedCPU()); nothing)
f()

using ProfileView, Profile, StaticArrays

Profile.clear()
@btime sim!($output, Ruleset());
@btime sim!($output, Ruleset(); proc=ThreadedCPU());
@btime sim!($output, $localdisp; proc=ThreadedCPU());
@btime sim!($output, $localdisp; proc=ThreadedCPU(), opt=SparseOpt());

@btime sim!($cu_output, $localdisp; proc=CuGPU());
@btime sim!($output, $growth);
@btime sim!($output, $localdisp, $allee);
@btime sim!($output, $growth; proc=CuGPU());
@btime sim!($output, $localdisp, $growth);
@btime sim!($output, Chain($localdisp, $growth));
@btime sim!($output, Chain($localdisp, $growth); proc=ThreadedCPU());
@btime sim!($cu_output, Chain($localdisp, $growth); proc=CuGPU());

@btime sim!($output, $wind, $allee, $growth);
@btime sim!($output, $wind, $allee, $growth; opt=SparseOpt());
@btime sim!($output, $wind, $allee, $growth; proc=ThreadedCPU());
@btime sim!($output, $wind, $allee, $growth; proc=ThreadedCPU(), opt=SparseOpt());

@btime sim!($output, Chain($localdisp, $growth), $wind);
@btime sim!($output, Chain($localdisp, $growth), $wind; proc=ThreadedCPU());

@btime sim!($output, ($localdisp, $growth, $allee));
@btime sim!($output, $localdisp, Chain($allee, $growth));
@btime sim!($output, $localdisp, Chain($growth, $allee));
@btime sim!($output, ($localdisp, $growth, $allee); opt=SparseOpt());
@btime sim!($output, Chain($localdisp, $growth, $allee); opt=SparseOpt());
@btime sim!($output, $localdisp, Chain($growth, $allee); opt=SparseOpt());
@btime sim!($output, ($localdisp, $growth, $allee); proc=ThreadedCPU());
@btime sim!($output, Chain($localdisp, $growth, $allee); proc=ThreadedCPU());
@btime sim!($output, ($localdisp, $growth, $allee); proc=ThreadedCPU(), opt=SparseOpt());
@btime sim!($output, $localdisp, Chain($growth, $allee); proc=ThreadedCPU(), opt=SparseOpt());
@btime sim!($cu_output, ($localdisp, $growth, $allee); proc=CuGPU());
@btime sim!($cu_output, Chain($localdisp, $growth, $allee); proc=CuGPU());

@btime sim!($output, $randomgrid, $gpu_wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism);
@btime sim!($output, $randomgrid, $gpu_wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism; opt=SparseOpt());
@btime sim!($output, $randomgrid, $gpu_wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism; proc=ThreadedCPU());
@btime sim!($output, $wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism; proc=ThreadedCPU());
@btime sim!($output, $wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism; opt=SparseOpt(), proc=ThreadedCPU());
@btime sim!($output, $randomgrid, $gpu_wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism; opt=SparseOpt(), proc=ThreadedCPU());
@btime sim!($output, $randomgrid, $gpu_wind, $localdisp, Chain($allee, $growth), $localdispP, Chain($alleeP, $parasitism); opt=SparseOpt(), proc=ThreadedCPU());
@btime sim!($output, $randomgrid, $gpu_wind, $localdisp, $allee, $growth, $localdispP, $alleeP, $parasitism; proc=CuGPU());
@btime sim!($cu_output, $randomgrid, $gpu_wind, $localdisp, $localdispP, Chain($allee, $growth, $alleeP, $parasitism); proc=CuGPU());
@btime sim!($output, Ruleset(); proc=CuGPU());
@btime sim!($output, Ruleset());

cu_output = adapt(CuArray, output);
sim!(cu_output, randomgrid, gpu_wind, localdisp, Chain(allee, growth), localdispP, alleeP, parasitism; proc=DynamicGrids.CuGPU());

using Adapt
@profview sim!(cu_output, Ruleset(); proc=CuGPU());
@btime sim!(cu_output, Ruleset(); proc=CuGPU());

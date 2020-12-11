using BenchmarkTools, Plots, KernelAbstractions, CUDA, Random

function setupsim(rules;
    tspan=DateTime(2020, 1, 1):Day(1):DateTime(2020, 4, 10),
    size=(100, 100), opt=NoOpt(), proc=SingleCPU(),
    growthmin=-5.0f0, growthmax=0.2f0
)
    # gr = rand(Float32, size...) .* (growthmax - growthmin) .- growthmin
    # init = (; H=rand(Bool, size) .* 100.0f0, rand=rand(Float32, size))
    hostpop .= 100.0f0
    rnge = Base.OneTo.(size) 
    o = ArrayOutput((; H=hostpop[rnge...], P=parapop[rnge...], rand=rand(Float32, size)); 
        aux=(; rH=r_host_aus[rnge..., :], rP=r_para_aus[rnge..., :]), 
        mask=Array(mask)[rnge...],
        tspan=tspan
    )
    rs = Ruleset(rules; opt=opt, proc=proc)
    # sd = DynamicGrids.SimData(DynamicGrids.extent(o), rs)
    # Warmup compilation
    # sim!(o, rs)
    return (; o, rs)
end
function bench(rules)
    suite = BenchmarkGroup()
    for opt in keys(opts)
        suite[opt] = BenchmarkGroup(["proc", "size"])
        for proc in keys(procs), s in sizes
            suite[opt][proc, s] = @benchmarkable sim!(x[:o], x[:rs]) setup=(x = setupsim($rules; size=($s, $s), proc=$(procs[proc]), opt=$(opts[opt])))
        end
    end
    results = run(suite, verbose=true)
end

procs = Dict("gpu"=>CuGPU(), "single"=>SingleCPU(), "threaded"=>ThreadedCPU())
opts = Dict("noopt"=>NoOpt(), "sparseopt"=>SparseOpt())
sizes = (100, 200, 300, 400)#, 1000, 2000)
rulegroups = (Wind=(randomgrid, gpu_wind, growth), Local=(localdisp, growth), Combined=(localdisp, randomgrid, gpu_wind, growth))
results = map(bench, rulegroups)
results[:Combined][@tagged 400]
minimum(results[:Combined].data["sparseopt"].data[("gpu", 400)].times) / 1e9

output = ArrayOutput((; H=parent(hostpop), P=parent(parapop), rand=rand(Float32, size(hostpop))); 
    aux=(; rH=r_host_aus, rP=r_para_aus), 
    mask=Array(mask),
    tspan=tspan
)

for key in keys(rulegroups), p in values(procs), o in values(opts)
    @show key p o
    @time sim!(output, rulegroups[key]; proc=p)
end

function plotbench(b, key)
    p = plot(; )
    for prockey in Tuple(keys(procs))
        if prockey != "gpu" # no sparseopt on GPU
            so_times = map(sizes) do s
                minimum(b.data["sparseopt"].data[(prockey, s)].times)  / 1e9
            end
            @show so_times
            plot!(p, collect(sizes), collect(so_times); label="Sparse $prockey", pallete=:Dark2_5, opacity=0.8)
        end
        no_times = map(sizes) do s
            minimum(b.data["noopt"].data[(prockey, s)].times) / 1e9
        end
        @show no_times
        plot!(p, collect(sizes), collect(no_times); label=uppercasefirst(prockey), pallete=:Dark2_5, opacity=0.8)
    end
    plot(p; 
        title=string(key), 
        ylims=(0, 0.4), 
        legend=:top, 
        tickfontcolor=RGB(0.2),
        legendfontcolor=RGB(0.2),
        minorgrid=false,
        xlabel="Grid side",
        ylabel="Time in seconds",
    )
end
pyplot()
theme(:vibrant)
using Plots: px
plot(map(plotbench, results, keys(results))...; layout=(1, 3), size=(1200, 400))

savefig("benchmarks.png")


### Other ########################################33

using ConstructionBase
ConstructionBase.constructorof(::Type{<:LinRange}) = _linrange
_linrange(start, stop, len, lendiv) = LinRange(start, stop, len)

hostpop_disag, parapop_disag, r_host_aus_disag, r_para_aus_disag, mask_disag = 
    map((hostpop, parapop, r_host_aus, r_para_aus, mask)) do A
        A = set(A, Lat=NoIndex(), Lon=NoIndex())
        # A = hasdim(A, Ti) ? set(A, Ti=NoIndex()) : A
        GeoData.disaggregate(Center(), A, (Lat(4), Lon(4), Ti(1)))
    end

A = set(r_host_aus, Lat=NoIndex(), Lon=NoIndex(), Ti=NoIndex())
typeof(index(r_host_aus_disag))

output = ArrayOutput((; H=hostpop_disag, P=parapop_disag, rand=rand(Float32, size(hostpop_disag))); 
    aux=(; rH=r_host_aus_disag, rP=r_para_aus_disag), mask=Array(mask_disag), tspan=tspan)

@time sim!(output, growth);
@time sim!(output, growth; proc=ThreadedCPU());
@time sim!(output, growth; proc=CuGPU());
@time sim!(output, localdisp);
@time sim!(output, localdisp; proc=ThreadedCPU());
@time sim!(output, localdisp; proc=CuGPU());
@time sim!(output, localdisp, growth);
@time sim!(output, localdisp, growth; proc=ThreadedCPU());
@time sim!(output, localdisp, growth; proc=CuGPU());
@time sim!(output, gpu_wind, randomgrid, growth);
@time sim!(output, gpu_wind, randomgrid, growth; proc=ThreadedCPU());
@time sim!(output, gpu_wind, randomgrid, growth; proc=CuGPU());
@time sim!(output, localdisp, gpu_wind, randomgrid, growth);
@time sim!(output, localdisp, gpu_wind, randomgrid, growth; proc=ThreadedCPU());
@time sim!(output, localdisp, gpu_wind, randomgrid, growth; proc=CuGPU());

@time sim!(output, growth, opt=SparseOpt());
@time sim!(output, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, growth; proc=CuGPU(), opt=SparseOpt());
@time sim!(output, localdisp, growth, opt=SparseOpt());
@time sim!(output, localdisp, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, localdisp, growth; proc=CuGPU(), opt=SparseOpt());
@time sim!(output, gpu_wind, randomgrid, growth, opt=SparseOpt());
@time sim!(output, gpu_wind, randomgrid, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, gpu_wind, randomgrid, growth; proc=CuGPU(), opt=SparseOpt());
@time sim!(output, localdisp, gpu_wind, randomgrid, growth, opt=SparseOpt());
@time sim!(output, localdisp, gpu_wind, randomgrid, growth; proc=ThreadedCPU(), opt=SparseOpt());
@time sim!(output, localdisp, gpu_wind, randomgrid, growth; proc=CuGPU(), opt=SparseOpt());

ruleset = (randomgrid, gpu_wind, localdisp, para_localdisp, host_para_growth)
@time sim!(output, ruleset);
@time sim!(output, ruleset; proc=ThreadedCPU());
@time sim!(output, ruleset; proc=CuGPU());
# using ProfileView, Profile
# Profile.clear()
# f(output, growth) = @profile for i in 1:100 sim!(output, randomgrid, gpu_wind, growth) end
# f(output, growth)
# ProfileView.view()

# output = ArrayOutput((; host=hostpop); tspan=tspan)
# error()

# @time sim!(output, localdisp);
# @btime sim!($output, $localdisp);
# @btime sim!($output, $localdisp; proc=ThreadedCPU());
# @btime sim!($output, $localdisp; proc=CuGPU());
# output = ArrayOutput((; host=hostpop); aux=(; growthrates=growthrates1), tspan=tspan)
# @btime sim!($output, Chain($localdisp, $growth));
# @btime sim!($output, Chain($localdisp, $growth); proc=ThreadedCPU());
# @btime sim!($output, Chain($localdisp, $growth); proc=CuGPU());
# # * use a `let` block for 0.2 so the variable is locally-scope

# output = ArrayOutput((; host=hostpop); aux=(; growthrates=growthrates1), tspan=tspan)
# sim!(output, (localdisp, wind, growth))
# @btime sim!($output, Chain($localdisp, $growth), $wind);
# @btime sim!($output, Chain($localdisp, $growth), $wind; proc=ThreadedCPU());
# @btime sim!($output, Chain($localdisp, $growth), $wind; proc=CuGPU());
# output = REPLOutput((; host=hostpop);tspan=tspan)
# output = ArrayOutput((; host=hostpop);tspan=tspan)
# sim!(output, wind)
# @btime sim!($output, $wind);
# @btime sim!($output, $wind; proc=ThreadedCPU());
# # THIS CRASHES @btime sim!($output, $wind; proc=CuGPU());

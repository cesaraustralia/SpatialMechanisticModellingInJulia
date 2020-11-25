using BenchmarkTools

#### Benchmarks ##########################################################

# Wind dispersal model For GPU randomisation

# As we don't have `rand()` inside the kernel, we have to calulate it 
# in another grid from outside a kernel, using a `Grid` rul.e
# We want rand to be padded with a neighborhood so we can sample
# extra values.

using Random
randomgrid = Grid{Tuple{},:rand}() do w
    if w isa CuArray
        CUDA.rand!(w)
    else
        Random.rand!(w)
    end
end
# init = zeros(Float32, 1000, 1000)
# o = ResultOutput((; rand=zero(hostpop)); tspan=1:100)
# @time sim!(o, randomgrid; proc=CuGPU())

# @btime sim!($o, $randomgrid; proc=CuGPU())
# @time sim!(o, randomgrid)
# @btime sim!($o, $randomgrid)

gpu_wind = SetNeighbors{Tuple{:rand,:H},Tuple{:rand,:H}}() do data, hood, index, rand, state
    # Ignore empty cells
    state > zero(state) || return nothing
    # Randomise a dispersal event using cell value of rand,
    # calculated in the rand Grid rule
    # rand < 0.2 || return nothing
    # Randomise a destination
    # spotrange = 10
    # Get random values for jump from the surrounding grid
    # randi = (index[1], index[2] + 1)
    # randj = (index[1] + 1, index[2] + 1)
    # jump = round.(Int, (data[:rand][randi...], data[:rand][randj...]) .* 2spotrange .- spotrange)
    # Make sure the destination is on the grid
    # jumpdest, is_inbounds = inbounds(jump .+ index, data[:host])
    # Update spotted cell if it's on the grid
    # is_inbounds && @inbounds add!(data[:host], state, jumpdest...)
    return nothing
end

# output = ArrayOutput((; host=hostpop, rand=zero(hostpop)); tspan=tspan)
# sim!(output, (randomgrid, gpu_wind))
# @btime sim!($output, $wind);
# @btime sim!($output, $wind; proc=ThreadedCPU());
# @btime sim!($output, $randomgrid, $gpu_wind; proc=CuGPU());


function setupsim(rules; 
    tspan=DateTime(2020, 1, 1):Day(1):DateTime(2020, 4, 10), 
    size=(100, 100), opt=NoOpt(), proc=SingleCPU(), 
    growthmin=-5.0, growthmax=0.2
)
    growthrates = rand(size...) .* (growthmax - growthmin) .- growthmin
    init = (; H=zeros(size...))
    o = ResultOutput(init; 
        aux=(; rH=growthrates),
        mask=parent(mask),
        tspan=tspan,
    )
    rs = Ruleset(rules; opt=opt, proc=proc)
    sd = DynamicGrids.SimData(DynamicGrids.extent(o), rs)
    # Warmup compilation
    sim!(o, rs; simdata=sd)
    return (; o, rs, sd)
end
function bench(rules)
    # procs = Dict("gpu"=>CuGPU(), "single"=>SingleCPU(), "threaded"=>ThreadedCPU())
    procs = Dict("single"=>SingleCPU(), "threaded"=>ThreadedCPU())
    opts = Dict("noopt"=>NoOpt(), "sparseopt"=>SparseOpt())
    sizes = (100, 200, 400)#, 1000, 2000)
    suite = BenchmarkGroup()
    suite[:noopt] = BenchmarkGroup()
    suite[:sparseopt] = BenchmarkGroup()
    for opt in keys(opts)
        suite[opt] = BenchmarkGroup(["proc", "size"])
        for proc in keys(procs), s in sizes
            suite[opt][proc, s] = @benchmarkable sim!(x[:o], x[:rs]; simdata=x[:sd]) setup=(x = setupsim($rules; size=($s, $s), proc=procs[$proc], opt=opts[$opt]))
        end
    end
    results = run(suite, verbose=true)
end
function plotbench(b)
    ps1 = map(keys(procs)) do prockey
        times = map(sizes) do s
            minimum(b.data["sparseopt"].data[(prockey, s)].times)
        end
        plot!(collect(sizes), collect(times) .* 10e-6; label="sparse $prockey")
    end
    ps2 = map(keys(procs)) do prockey
        times = map(sizes) do s
            minimum(b.data["noopt"].data[(prockey, s)].times)
        end
        plot!(collect(sizes), collect(times) .* 10e-6; label=prockey)
    end
    plot(ps1..., ps2...; yaxis=:log);
end
rulegroups = (wind=(wind, growth), loc=(localdisp, growth), combined=(localdisp, wind, growth))
results = map(bench, rulegroups)
# results[@tagged 500]
# plotbench(results)



### Other ########################################33


# output = ArrayOutput((; host=hostpop); aux=(; rH=growthrates1), tspan=tspan)
# sim!(output, (localdisp, growth, viewer))
# sim!(output, (wind, growth, viewer))

# @time sim!(output, growth);
# @time sim!(output, localdisp, growth);
# @time sim!(output, wind, localdisp, growth);
# @time sim!(output, wind, growth);
# @time sim!(output, localdisp, growth; proc=ThreadedCPU());
# @time sim!(output, wind, growth; proc=ThreadedCPU());
# @time sim!(output, localdisp, growth; proc=CuGPU());
# @time sim!(output, wind, growth; proc=CuGPU());

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

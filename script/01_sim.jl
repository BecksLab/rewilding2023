println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

first_sim = parse(Int, ARGS[1])
last_sim = parse(Int, ARGS[2])

println("Running parameters combination from $(ARGS[1]) to $(ARGS[2])")

#Flag enables all the workers to start on the project of the current dir
#dir = pwd() * "/"
dir = expanduser("~/rewilding2023/")
flag = "--project=" * dir
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(10 - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -2) cores")

@everywhere import Pkg
@everywhere using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
@everywhere include("../src/foodweb_building.jl")

@everywhere using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
@everywhere using DifferentialEquations, SparseArrays
@everywhere using EcologicalNetworksDynamics
@everywhere using Random, Plots, Distributions, DataFrames
@everywhere using ProgressMeter
@everywhere include("../src/misc.jl")
@everywhere include("../src/simulation_methods.jl")
@everywhere include("../src/foodweb_measure.jl")

import Random.seed!

seed!(22)
###################
#  CS experiment  #
###################

# Idea: same output for the three experiments
# Big df result with a "step" key: "present", "extirpated", "reintroduced"

# Load Foodweb
fw_comb_df = DataFrame(Arrow.Table(joinpath(dir, "data/sim_param.arrow")))

if last_sim > size(fw_comb_df, 1)
    last_sim = size(fw_comb_df, 1)
end
println("Running param sim from lines $first_sim to $last_sim")

# Reshape arrays and make a vector
fw_comb_df[!, :fw] = map(x -> reshape_array(x), fw_comb_df[:, :fw])
fw_comb = NamedTuple.(eachrow(fw_comb_df))
#fw_comb = fw_comb[repeat([27], 50)]

# How many timestep to keep:
last_timestep = 500
burn_in_timestep = 2000
extinct_threshold = 1e-5
ndigits_bm = 5
gc_prob = .01

sim = @showprogress pmap(x -> sim_pred_present_extirpation_reintroduction(x;
                                                                          species_names = nothing,
                                                                          last_timestep = last_timestep,
                                                                          burn_in_timestep = burn_in_timestep,
                                                                          extinct_threshold = extinct_threshold,
                                                                          ndigits_bm = n_digits,
                                                                          gc_prob = gc_prob
                                                                         ),
                         fw_comb[first_sim:last_sim]; on_error = ex -> missing,
                         batch_size = 100
                        )

sim_df = reduce(vcat, skipmissing(sim))
println("Finished sim")
file = string("data/sim_", first_sim, "_", last_sim, ".arrow")
Arrow.write(joinpath(dir, file), sim_df)

println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

#Flag enables all the workers to start on the project of the current dir
#dir = pwd() * "/"
dir = expanduser("~/rewilding2023/")
flag = "--project=" * dir
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(20 - 1, exeflags=flag)
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

function stoch_gc(gc_thre)
    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
        GC.safepoint()
    end
end

sim = @showprogress pmap(x -> merge(
                                    (sim_id = x.sim_id, fw = x.fw, Z = x.Z, nb_link = x.nb_link),
                                    (
                                     out = begin
                                         fw = FoodWeb(x.fw, Z = x.Z)
                                         p = ModelParameters(fw)
                                         B0 = Base.rand(richness(fw))
                                         m = sim_steady_state_last(p, B0,
                                                                   last = last_timestep,
                                                                   burn_in = burn_in_timestep,
                                                                   extinction_threshold = extinct_threshold
                                                                  )
                                         scenario_output(m, B0, fw;
                                                         last = last_timestep,
                                                         scenario = "pred_present",
                                                         i_extirpated = missing,
                                                         tlvl_extirpated = missing,
                                                         i_introduced = missing,
                                                         tlvl_introduced = missing,
                                                         n_digits = ndigits_bm
                                                        )
                                     end,
                                    ).out
                                   ),
                         fw_comb; on_error = ex -> missing,
                         batch_size = 100
                        )
sim_df = DataFrame(skipmissing(sim))
Arrow.write(joinpath(dir, "data/sim_pred_present.arrow"), sim_df)

# Add check for dead species: put their biomass to 0
# Add check for negative bm: put their biomass to 0
#
# Without top predator
sim_extinction = @showprogress pmap(x -> merge(
                                               (sim_id = x.sim_id, fw = x.fw, Z = x.Z, nb_link = x.nb_link),
                                               (out = begin
                                                    stoch_gc(gc_prob)
                                      fw = FoodWeb(x.fw, Z = x.Z)
                                      p = ModelParameters(fw)
                                      # Max tlvl among alive species
                                      max_tlvl_alive = findmax(x.tlvl[x.species_alive])
                                      extirpated_tl = max_tlvl_alive[1]
                                      # species index in the original community
                                      extirpated_i = x.species_alive[max_tlvl_alive[2]]
                                      init_bm = sanatize_biomass(x.bm_sp, x.species_alive)
                                      init_bm[introduced_i] = 0
                                      m = sim_steady_state_last(p, init_bm,
                                                                last = last_timestep,
                                                                burn_in = burn_in_timestep,
                                                                extinction_threshold = extinct_threshold
                                                               )
                                      scenario_output(m, init_bm, fw;
                                                      last = last_timestep,
                                                      scenario = "pred_extirpated",
                                                      i_extirpated = extirpated_i,
                                                      tlvl_extirpated = extirpated_tl,
                                                      i_introduced = missing,
                                                      tlvl_introduced = missing,
                                                      n_digits = ndigits_bm
                                                     )
                                  end,
                                 ).out
                                ),
                                    sim; on_error = ex -> missing,
                                    batch_size = 100
                                   )

sim_extinction_df = DataFrame(skipmissing(sim_extinction))
Arrow.write(joinpath(dir, "data/sim_extinction.arrow"), sim_extinction_df)

# Re-introduction:
sim_reintroduction = pmap(x -> merge(
                                     (sim_id = x.sim_id, fw = x.fw, Z = x.Z, nb_link = x.nb_link),
                                     (out = begin
                                          stoch_gc(gc_prob)
                                          fw = FoodWeb(x.fw, Z = x.Z)
                                          # Remove or add link to the predator
                                          if (!ismissing(x.nb_link))
                                              f = add_remove_link_species(fw,
                                                                          x.i_extirpated;
                                                                          nb_link = x.nb_link)
                                          else
                                              f = fw.A
                                          end
                                          fw = FoodWeb(f, Z = x.Z)
                                          p = ModelParameters(fw)
                                          init_bm = sanatize_biomass(x.bm_sp, x.species_alive)
                                          # Re-introduce the predator:
                                          init_bm[x.i_extirpated] = mean(init_bm[init_bm .> 0])
                                          m = sim_steady_state_last(p, init_bm,
                                                                last = last_timestep,
                                                                burn_in = burn_in_timestep,
                                                                extinction_threshold = extinct_threshold
                                                                )
                                          scenario_output(m, init_bm, fw;
                                                          last = last_timestep,
                                                          scenario = "pred_reintroduced",
                                                          i_extirpated = missing,
                                                          tlvl_extirpated = missing,
                                                          i_introduced = x.i_extirpated,
                                                          tlvl_introduced = x.tlvl_extirpated,
                                                          n_digits = ndigits_bm
                                           )
                                      end,
                                     ).out
                                    ),
                          sim_extinction; on_error = ex -> missing,
                          batch_size = 100
                         )
sim_reintroduction_df = DataFrame(skipmissing(sim_reintroduction))
Arrow.write(joinpath(dir, "data/sim_reintroduction.arrow"), sim_reintroduction_df)

sim_tot = [sim; sim_extinction; sim_reintroduction]
sim_tot_df = DataFrame(skipmissing(sim_tot))

Arrow.write(joinpath(dir, "data/sim_total.arrow"), sim_tot_df)

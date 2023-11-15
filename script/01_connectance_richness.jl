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
addprocs(15 - 1, exeflags=flag)
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

sim = @showprogress pmap(x -> merge(
                                    (sim_id = x.sim_id, fw = x.fw, Z = x.Z, nb_link = x.nb_link),
                                    (
                                     out = begin
                                         fw = FoodWeb(x.fw, Z = x.Z)
                                         p = ModelParameters(fw)
                                         B0 = Base.rand(richness(fw))
                                         m = sim_steady_state_last(p, B0, last = 100)
                                         scenario_output(m, B0, fw;
                                                         last = 100,
                                                         scenario = "pred_present",
                                                         i_extirpated = missing,
                                                         tlvl_extirpated = missing,
                                                         i_introduced = missing,
                                                         tlvl_introduced = missing,
                                                        )
                                     end,
                                    ).out
                                   ),
                         fw_comb; on_error = ex -> missing,
                         batch_size = 100
                        )
sim_df = DataFrame(skipmissing(sim))
Arrow.write(joinpath(dir, "data/sim_pred_present.arrow"), sim_df)


# Without top predator
sim_extinction = @showprogress pmap(x -> merge(
                                               (sim_id = x.sim_id, fw = x.fw, Z = x.Z, nb_link = x.nb_link),
                                               (out = begin
                                      fw = FoodWeb(x.fw, Z = x.Z)
                                      p = ModelParameters(fw)
                                      # Max tlvl among alive species
                                      max_tlvl_alive = findmax(x.tlvl[x.species_alive])
                                      introduced_tl = max_tlvl_alive[1]
                                      # species index in the original community
                                      introduced_i = x.species_alive[max_tlvl_alive[2]]
                                      init_bm = x.bm_sp
                                      init_bm[introduced_i] = 0
                                      m = sim_steady_state_last(p, init_bm, last = 100)
                                      scenario_output(m, init_bm, fw;
                                                      last = 100,
                                                      scenario = "pred_extirpated",
                                                      i_extirpated = introduced_i,
                                                      tlvl_extirpated = introduced_tl,
                                                      i_introduced = missing,
                                                      tlvl_introduced = missing,
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
                                          init_bm = x.bm_sp
                                          # Re-introduce the predator:
                                          init_bm[x.i_extirpated] = mean(init_bm[init_bm .> 0])
                                          m = sim_steady_state_last(p, init_bm, last = 100)
                                          scenario_output(m, init_bm, fw;
                                                          last = 100,
                                                          scenario = "pred_reintroduced",
                                                          i_extirpated = missing,
                                                          tlvl_extirpated = missing,
                                                          i_introduced = x.i_extirpated,
                                                          tlvl_introduced = x.tlvl_extirpated
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

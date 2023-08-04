
println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

#Flag enables all the workers to start on the project of the current dir
dir = "~/rewilding2023/"
flag = "--project=" * dir
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -2) cores")

@everywhere import Pkg
@everywhere using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
@everywhere include(dir * "src/foodweb_building.jl")

@everywhere using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
@everywhere using DifferentialEquations, SparseArrays
@everywhere using EcologicalNetworksDynamics
@everywhere using Random, Plots, Distributions, DataFrames
@everywhere include(dir * "src/misc.jl")
@everywhere include(dir * "src/simulation_methods.jl")
@everywhere include(dir * "src/foodweb_measure.jl")

import Random.seed!

seed!(22)
###################
#  CS experiment  #
###################

# Idea: same output for the three experiments
# Big df result with a "step" key: "present", "extirpated", "reintroduced"

# Load Foodweb
fw_comb_df = DataFrame(Arrow.Table(dir * "data/fw_C_S.arrow"))

# Reshape arrays and make a vector
fw_comb_df[!, :fw] = map(x -> reshape_array(x), fw_comb_df[:, :fw])
fw_comb = NamedTuple.(eachrow(fw_comb_df))


sim = pmap(x -> merge(
                      (fw_id = x.fw_id,),
                      (
                       out = begin
                           fw = FoodWeb(x.fw, Z = 100)
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
           fw_comb
          )

sim_df = DataFrame(sim)

# Without top predator
sim_extinction = pmap(x -> merge(
                                 (fw_id = x.fw_id, ),
                                 (out = begin
                                      fw = FoodWeb(x.A_alive, Z = 100)
                                      p = ModelParameters(fw)
                                      introduced_tl, introduced_i = findmax(x.tlvl)
                                      init_bm = x.bm_sp[x.species_alive]
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
                      sim
                     )


# Re-introduction:
sim_reintroduction = pmap(x -> merge(
                                     (fw_id = x.fw_id, ),
                                     (out = begin
                                          fw = FoodWeb(x.A_init, Z = 100)
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
                          sim_extinction
                         )

sim_tot = [sim; sim_extinction; sim_reintroduction]
sim_tot_df = DataFrame(sim_tot)

Arrow.write(dir * "data/simCS.arrow",sim_tot_df)

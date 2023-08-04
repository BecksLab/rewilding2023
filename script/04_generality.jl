using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
using DifferentialEquations, SparseArrays

using EcologicalNetworksDynamics
using Random, Plots, Distributions, DataFrames
using Distributed
include("src/misc.jl")
include("src/simulation_methods.jl")
include("src/foodweb_measure.jl")

###################
#  CS experiment  #
###################

# Idea: same output for the three experiments
# Big df result with a "step" key: "present", "extirpated", "reintroduced"

# Load Foodweb
fw_comb_df = DataFrame(Arrow.Table("data/fw_G.arrow"))

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
                                          # Remove or add link to the predator
                                          f = add_remove_link_species(fw,
                                                                      x.i_extirpated;
                                                                      nb_link = x.nb_link)
                                          fw = FoodWeb(f, Z = 100)
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

Arrow.write("data/simG.arrow",sim_tot_df)

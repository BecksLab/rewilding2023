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
fw_comb_df = DataFrame(Arrow.Table("data/fw_C_S.arrow"))

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
                           troph_class_init = trophic_classes(x.fw)
                           m = sim_steady_state_last(p, B0, last = 100)
                           merge(
                                 sim_output(m, last = 100),
                                 (
                                  richness_init = richness(fw),
                                  bm_init = B0,
                                  bm_init_pred = sum(B0[troph_class_init.top_predators]),
                                  bm_init_cons = sum(B0[troph_class_init.intermediate_consumers]),
                                  bm_init_prod = sum(B0[troph_class_init.producers]),
                                 )
                                )
                       end,
                      ).out
                     ),
           fw_comb[1:2]
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
                                      troph_class_init = trophic_classes(x.A_alive)
                                      m = sim_steady_state_last(p, init_bm, last = 100)
                                      merge(
                                            sim_output(m, last = 100),
                                            (
                                             richness_init = richness(fw),
                                             bm_init = init_bm,
                                             bm_init_pred = sum(init_bm[troph_class_init.top_predators]),
                                             bm_init_cons = sum(init_bm[troph_class_init.intermediate_consumers]),
                                             bm_init_prod = sum(init_bm[troph_class_init.producers]),
                                             tlvl_extinct = introduced_tl,
                                             i_extinct = introduced_i,
                                             A_init = x.A_alive
                                            )
                                           )
                                  end,
                                 ).out
                                ),
                      sim[1:2]
                     )


# Re-introduction:
sim_reintroduction = pmap(x -> merge(
                                     (fw_id = x.fw_id, ),
                                     (out = begin
                                          fw = FoodWeb(x.A_init, Z = 100)
                                          p = ModelParameters(fw)
                                          init_bm = x.bm_sp
                                          # Re-introduce the predator:
                                          init_bm[x.i_extinct] = mean(init_bm[init_bm .> 0])
                                          troph_class_init = trophic_classes(x.A_init)
                                          m = sim_steady_state_last(p, init_bm, last = 100)
                                          merge(
                                                sim_output(m, last = 100),
                                                (
                                                 richness_init = x.richness,
                                                 bm_init = init_bm,
                                                 bm_init_pred = sum(init_bm[troph_class_init.top_predators]),
                                                 bm_init_cons = sum(init_bm[troph_class_init.intermediate_consumers]),
                                                 bm_init_prod = sum(init_bm[troph_class_init.producers]),
                                                )
                                               )
                                      end,
                                     ).out
                                    ),
                          sim_extinction[1:2]
                         )

# Z experiment
Z = [10, 20, 50, 100, 200, 500, 1000]

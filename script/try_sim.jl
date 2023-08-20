
fw_comb_test = NamedTuple.(eachrow(filter(:rep => x -> x == 1, fw_comb_df)))
sim = []
for x in fw_comb_test
    println("$(x.fw_id)")
    fw = FoodWeb(x.fw, Z = 100)
    p = ModelParameters(fw)
    B0 = Base.rand(richness(fw))
    m = sim_steady_state_last(p, B0, last = 100)
    out = scenario_output(m, B0, fw;
                          last = 100,
                          scenario = "pred_present",
                          i_extirpated = missing,
                          tlvl_extirpated = missing,
                          i_introduced = missing,
                          tlvl_introduced = missing,
                         )
    push!(sim, merge((fw_id = x.fw_id,), out))
end

sim_extinction = []
x = sim[2]
for x in sim
    println("$(x.fw_id)")
    fw = FoodWeb(x.A_alive, Z = 100)
    p = ModelParameters(fw)
    introduced_tl, introduced_i = findmax(x.tlvl)
    init_bm = x.bm_sp[x.species_alive]
    init_bm[introduced_i] = 0
    m = sim_steady_state_last(p, init_bm, last = 100)
    out = scenario_output(m, init_bm, fw;
                          last = 100,
                          scenario = "pred_extirpated",
                          i_extirpated = introduced_i,
                          tlvl_extirpated = introduced_tl,
                          i_introduced = missing,
                          tlvl_introduced = missing,
                         )
    push!(sim_extinction, merge((fw_id = x.fw_id,), out))
end

steady_sim = simulate(p, init_bm,
                      callback = CallbackSet(
                                             TerminateSteadyState(1e-6, 1e-4),
                                             ExtinctionCallback(1e-6, p, true)
                                            )
                     )
    last_biomass = steady_sim[:, end]

m2 =  simulate(p, last_biomass,
          tmax = 200,
          saveat = 0:1:200,
          callback = CallbackSet(
                                 ExtinctionCallback(1e-6, p, true)
                                ),
            )
m2[:, end]
sim_extinction[2].bm_sp

sim_reintroduction = []
for x in sim_extinction
    println("$(x.fw_id)")
    fw = FoodWeb(x.A_init, Z = 100)
    p = ModelParameters(fw)
    init_bm = x.bm_sp
    # Re-introduce the predator:
    init_bm[x.i_extirpated] = mean(init_bm[init_bm .> 0])
    m = sim_steady_state_last(p, init_bm, last = 100)
    out = scenario_output(m, init_bm, fw;
                    last = 100,
                    scenario = "pred_reintroduced",
                    i_extirpated = missing,
                    tlvl_extirpated = missing,
                    i_introduced = x.i_extirpated,
                    tlvl_introduced = x.tlvl_extirpated
                   )
    push!(sim_reintroduction, merge((fw_id = x.fw_id,), out))
end

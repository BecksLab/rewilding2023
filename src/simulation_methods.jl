#' Simulate until steady state and then simulate for a `last` number of timestep with a step
# of 1.
function sim_steady_state_last(p, B0;
        last = 100,
        burn_in = 100,
        extinction_threshold = 1e-6,
        verbose = true,
        kwargs...)
    steady_sim = simulate(p, B0, callback = CallbackSet(
                                                     TerminateSteadyState(1e-6, 1e-4),
                                                     ExtinctionCallback(extinction_threshold, p, verbose)
                                                    ),
                       kwargs...
                      )
    last_biomass = steady_sim[:, end]
    run_time = burn_in + last

    simulate(p, last_biomass,
             tmax = run_time,
             saveat = 0:1:run_time,
             callback = CallbackSet(
                                    ExtinctionCallback(extinction_threshold, p, verbose)
                                   ),
             kwargs...
            )
end

# The output extracted from the simulation
function sim_output(m; last = 100)
    p = get_parameters(m)
    bm = biomass(m, last = last)
    troph = trophic_structure(m, last = last)
    troph_class = trophic_classes(troph.alive_A)
    pref_alive = p.functional_response.ω[troph.alive_species, troph.alive_species]
    omni = omnivory(pref_alive, weighted = true)
    cv = coefficient_of_variation(m, last = last)
    int = empirical_interaction_strength(m, p, last = last).mean
    int_per_capita = int ./ bm.species

    (
     richness = richness(m),
     persistence = species_persistence(m),
     bm_sp = bm.species,
     bm_total = bm.total,
     bm_pred = sum(bm.species[troph_class.top_predators]),
     bm_cons = sum(bm.species[troph_class.intermediate_consumers]),
     bm_prod = sum(bm.species[troph_class.producers]),
     species_alive = troph.alive_species,
     A_alive = troph.alive_A,
     tlvl = troph.alive_trophic_level,
     tlvl_max = troph.max,
     tlvl_mean = troph.mean,
     tlvl_w_mean = troph.weighted_mean,
     omnivory = omni,
     omnivory_mean = mean(omni),
     cv_com = cv.community,
     cv_species = cv.species,
     cv_species_mean = cv.species_mean,
     int_mean = mean(int[int .> 0]),
     int_per_cap_mean = mean(int_per_capita[int_per_capita .> 0]),
    )
end

function scenario_output(m, B0, fw;
        last = 100,
        scenario = "pred_present",
        i_extirpated = missing,
        tlvl_extirpated = missing,
        i_introduced = missing,
        tlvl_introduced = missing
    )

    troph_class_init = trophic_classes(fw)
    merge(
          (scenario = scenario, ),
          sim_output(m, last = last),
          (
           richness_init = richness(B0),
           bm_init = B0,
           bm_init_pred = sum(B0[troph_class_init.top_predators]),
           bm_init_cons = sum(B0[troph_class_init.intermediate_consumers]),
           bm_init_prod = sum(B0[troph_class_init.producers]),
           tlvl_extirpated = tlvl_extirpated,
           i_extirpated = i_extirpated,
           i_introduced = i_introduced,
           tlvl_introduced = tlvl_introduced,
           A_init = fw.A
          )
   )
end

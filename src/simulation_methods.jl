#' Simulate until steady state and then simulate for a `last` number of timestep with a step
# of 1.
function sim_steady_state_last(p, B0;
        last = 500,
        burn_in = 1500,
        extinction_threshold = 1e-5,
        verbose = true,
        min_t_steady_state = 1000,
        gc_thre = .01,
        kwargs...)
    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
        GC.safepoint()
    end
    # steady_sim = simulate(p, B0, callback = CallbackSet(
                                                        # TerminateSteadyState(min_t = min_t_steady_state),
                                                        # ExtinctionCallback(extinction_threshold, p, verbose)
                                                       # ),
                       # kwargs...
                      # )
    #steady_sim
    #last_biomass = steady_sim[:, end]
    run_time = burn_in + last

    simulate(p, B0,
             tmax = run_time,
             saveat = 0:1:run_time,
             callback = CallbackSet(
                                    ExtinctionCallback(extinction_threshold, p, verbose)
                                   ),
             kwargs...
            )
end

# The output extracted from the simulation
function sim_output(m; last = 100, n_digits = 7)
    p = get_parameters(m)
    bm = biomass(m, last = last)
    p = get_parameters(m)
    fw = p.network.A
    troph = trophic_structure(m, last = last)
    troph_class = trophic_classes(fw)
    pref_alive = p.functional_response.ω[troph.alive_species, troph.alive_species]
    omni = omnivory(pref_alive, weighted = true)
    cv = coefficient_of_variation(m, last = last)
    int = empirical_interaction_strength(m, p, last = last).mean
    int_per_capita = int ./ bm.species
    A_alive = troph.alive_A
    (
     richness = richness(m),
     persistence = species_persistence(m),
     bm_sp = round.(bm.species, digits = n_digits),
     bm_total = bm.total,
     bm_pred = sum(bm.species[troph_class.top_predators]),
     bm_cons = sum(bm.species[troph_class.intermediate_consumers]),
     bm_prod = sum(bm.species[troph_class.producers]),
     species_alive = troph.alive_species,
     tlvl = trophic_levels(fw),
     tlvl_max = troph.max,
     tlvl_mean = troph.mean,
     tlvl_w_mean = troph.weighted_mean,
     connectance = connectance(A_alive),
     omnivory = omni,
     omnivory_mean = mean(omni),
     stab_com = 1 / cv.community,
     stab_species = 1 ./ cv.species,
     stab_species_mean = 1 / cv.species_mean,
    )
end

function scenario_output(m, B0, fw;
        last = 100,
        scenario = "pred_present",
        i_extirpated = missing,
        tlvl_extirpated = missing,
        i_introduced = missing,
        tlvl_introduced = missing,
        kwargs...
    )

    troph_class_init = trophic_classes(fw)
    merge(
          (scenario = scenario, ),
          sim_output(m; last = last, kwargs...),
          (
           richness_init = richness(B0),
           bm_init = B0,
           bm_init_pred = sum(B0[troph_class_init.top_predators]),
           bm_init_cons = sum(B0[troph_class_init.intermediate_consumers]),
           bm_init_prod = sum(B0[troph_class_init.producers]),
           tlvl_extirpated = tlvl_extirpated,
           i_extirpated = i_extirpated,
           i_introduced = i_introduced,
           tlvl_introduced = tlvl_introduced
          )
   )
end

"""

# Examples
sanatize_biomass([1, -1, 1], [1, 2, 3]) == [1, 0, 1]
sanatize_biomass([1, -1, 1], [1, 2]) == [1, 0, 0]
"""
function sanatize_biomass(bm, alive)
    # For negative biomass
    bm[bm .< 0] .= 0

    # For dead species
    in_alive = in(alive)
    mask_no_alive = .! in_alive.(1:length(bm))
    bm[mask_no_alive] .= 0
    bm
end

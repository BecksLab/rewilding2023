function test_correct_species_bm(sim_result)
    species_bm_pred_present = Dict(zip(sim_result.species[1], sim_result.bm_sp[1]))
    species_bm_init_extirpation = Dict(zip(sim_result.species[2], sim_result.bm_init[2]))
    species_bm_extirpation = Dict(zip(sim_result.species[2], sim_result.bm_sp[2]))
    species_bm_init_reintroduction = Dict(zip(sim_result.species[3], sim_result.bm_init[3]))

    compare_extirpation_pred_present = filter(
                                              p ->
                                              p.first !== sim_result.species_extirpated[2] &&
                                              any(p.first .== sim_result.species_alive[1]),
                                              species_bm_pred_present
                                             )
    test1 = compare_extirpation_pred_present == species_bm_init_extirpation

    compare_reintroduction_extirpation = filter(p -> p.first !== sim_result.species_introduced[3],
                                                species_bm_init_reintroduction)
    filtered_species_bm_extirpation = filter(p -> any(p.first .== sim_result.species_alive[2]),
                                             species_bm_extirpation)
    test2 = compare_reintroduction_extirpation == filtered_species_bm_extirpation

    test1 && test2
end

"""

# Test
## Second species is the top predator
ti = sim_pred_present_extirpation_reintroduction((sim_id = 1, fw = [0 0; 1 0], Z = 10, nb_link = missing,))
ti[2].species_extirpated == ti[3].species_introduced == "s2"
ti[1].richness == ti[3].richness == ti[2].richness + 1 == 2

## With changing species name
ti = sim_pred_present_extirpation_reintroduction((sim_id = 1, fw = [0 0; 1 0], Z = 10, nb_link = missing,),
                                                 species_names = ["carrot", "rabbit"]
                                                )
ti[2].species_extirpated == ti[3].species_introduced == "rabbit"
ti[1].richness == ti[3].richness == ti[2].richness + 1 == 2

## With with a three species foodweb
fw = [0 0 0; 1 0 0; 0 1 0]
ti = sim_pred_present_extirpation_reintroduction((sim_id = 1, fw = fw, Z = 10, nb_link = missing,),
                                                 species_names = ["carrot", "rabbit", "fox"]
                                                )
ti[2].species_extirpated == ti[3].species_introduced == "fox"
ti[1].richness == ti[3].richness == ti[2].richness + 1 == 3

## With with a four species foodweb
fw = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0]
ti = sim_pred_present_extirpation_reintroduction((sim_id = 1, fw = fw, Z = 10, nb_link = missing,),
                                                 species_names = ["carrot", "rabbit", "fox", "wolf"]
                                                )
ti[2].species_extirpated == ti[3].species_introduced == "wolf"
ti[1].richness == ti[3].richness == ti[2].richness + 1 == 4

## With with an extinction at warm-up
fw = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0]
ti = sim_pred_present_extirpation_reintroduction((sim_id = 1, fw = fw, Z = 1, nb_link = missing,),
                                                 species_names = ["carrot", "rabbit", "fox", "wolf"]
                                                )
ti[2].species_extirpated == ti[3].species_introduced == "fox"
ti[1].richness == ti[3].richness == ti[2].richness + 1 == 3

sim = @showprogress pmap(x -> sim_pred_present_extirpation_reintroduction(x),
                         [fw_comb[[930; 9000]]; (missing, )]; batch_size = 100, on_error = ex -> missing
                        )
reduce(vcat, skipmissing(sim))
"""
function sim_pred_present_extirpation_reintroduction(
        sim_data;
        species_names = nothing,
        last_timestep = 500,
        burn_in_timestep = 1500,
        extinct_threshold = 1e-5,
        ndigits_bm = 5,
        gc_prob = .01
    )

    #
    if isnothing(species_names)
        fw = FoodWeb(sim_data.fw, Z = sim_data.Z)
    else
        fw = FoodWeb(sim_data.fw, Z = sim_data.Z, species = species_names)
    end
    p = ModelParameters(fw)
    B0 = Base.rand(richness(fw))
    m = sim_steady_state_last(p, B0,
                              last = last_timestep,
                              burn_in = burn_in_timestep,
                              extinction_threshold = extinct_threshold,
                              gc_thre = gc_prob
                         )

    pred_present = scenario_output(m, B0, fw;
                                   last = last_timestep,
                                   scenario = "pred_present",
                                   species_extirpated = missing,
                                   tlvl_extirpated = missing,
                                   species_introduced = missing,
                                   tlvl_introduced = missing,
                                   n_digits = ndigits_bm
                                  )

    extirpated = sim_predator_extirpation(p, pred_present.species_alive,
                                          pred_present.bm_sp, pred_present.tlvl;
                                          last_timestep = last_timestep,
                                          burn_in_timestep = burn_in_timestep,
                                          extinct_threshold = extinct_threshold,
                                          gc_prob = gc_prob
                                         )
    reintroduction = sim_predator_reintroduction(p,
                                                 extirpated.species,
                                                 extirpated.species_alive,
                                                 extirpated.bm_sp,
                                                 extirpated.species_extirpated,
                                                 sim_data.nb_link;
                                                 last_timestep = last_timestep,
                                                 burn_in_timestep = burn_in_timestep,
                                                 extinct_threshold = extinct_threshold,
                                                 gc_prob = gc_prob
                                                )
    out = [merge((sim_id = sim_data.sim_id, ), i) for i in [pred_present; extirpated; reintroduction]]
    DataFrame(out)
end

function sim_predator_reintroduction(p,
        species,
        alive_species,
        bm_species,
        extirpated_species,
        nb_link;
        last_timestep = 500,
        burn_in_timestep = 1500,
        extinct_threshold = 1e-5,
        ndigits_bm = 5,
        gc_prob = .01
    )
    # Retrieve species name and position
    total_species = p.network.species
    species_idxs = Dict(zip(total_species, 1:length(total_species)))
    alive = alive_species
    # Species to keep (alive species + extirpated predator)
    species_to_keep = [alive; extirpated_species]
    idxs_to_keep = [species_idxs[i] for i in species_to_keep]

    # Remove or add link to the extirpated predator
    if (!ismissing(nb_link))
        f = add_remove_link_species(p.network,
                                    species_idxs[extirpated_species];
                                    nb_link = nb_link)
    else
        f = p.network.A
    end
    fw = FoodWeb(f[idxs_to_keep, idxs_to_keep],
                 M = p.network.M[idxs_to_keep],
                 species = species_to_keep,
                 metabolic_class = p.network.metabolic_class[idxs_to_keep]
                )
    p = ModelParameters(fw)
    # Begin with predator at the average biomass of all other species
    bm_species_dict = Dict(zip(species, bm_species))
    reordered_bm_species = [bm_species_dict[i] for i in alive]
    init_bm = [reordered_bm_species; mean(values(reordered_bm_species))]
    #@infiltrate
    m = sim_steady_state_last(p, init_bm,
                              last = last_timestep,
                              burn_in = burn_in_timestep,
                              extinction_threshold = extinct_threshold,
                              gc_thre = gc_prob
                             )
    scenario_output(m, init_bm, fw;
                    last = last_timestep,
                    scenario = "pred_reintroduced",
                    species_extirpated = missing,
                    tlvl_extirpated = missing,
                    species_introduced = extirpated_species,
                    tlvl_introduced = missing,
                    n_digits = ndigits_bm
                   )
end




# Extinction
function sim_predator_extirpation(p, alive_species, bm_species, tlvl;
        last_timestep = 500,
        burn_in_timestep = 1500,
        extinct_threshold = 1e-5,
        ndigits_bm = 5,
        gc_prob = .01
    )
    # Retrieve species name and position
    species = p.network.species
    species_idxs =  Dict(zip(species, 1:length(species)))
    alive = alive_species
    alive_idxs = [species_idxs[i] for i in alive]

    # Find the predator to extirpate
    max_tlvl_alive = findmax(tlvl[alive_idxs])
    # Species to extirpate
    extirpated_tl = max_tlvl_alive[1]
    extirpated_species = alive[max_tlvl_alive[2]]

    # Species to keep for the simulation (alive minus extirpated):
    species_to_keep = alive[alive .!== extirpated_species]
    idxs_to_keep = [species_idxs[i] for i in species_to_keep]

    # Simulation preparation
    ## Build the new FW
    fw = FoodWeb(p.network.A[idxs_to_keep, idxs_to_keep],
                  M = p.network.M[idxs_to_keep],
                  species = species_to_keep,
                  metabolic_class = p.network.metabolic_class[idxs_to_keep]
                 )
    ## Set parameters and initial biomass
    p = ModelParameters(fw)
    init_bm = bm_species[idxs_to_keep]
    ## Simulate
    m = sim_steady_state_last(p, init_bm,
                              last = last_timestep,
                              burn_in = burn_in_timestep,
                              extinction_threshold = extinct_threshold,
                              gc_thre = gc_prob
                             )
    # Compute output
    scenario_output(m, init_bm, fw;
                    last = last_timestep,
                    scenario = "pred_extirpated",
                    species_extirpated = extirpated_species,
                    tlvl_extirpated = extirpated_tl,
                    species_introduced = missing,
                    tlvl_introduced = missing,
                    n_digits = ndigits_bm
                   )

end

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
    species = p.network.species
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
     fw = fw,
     species = species,
     M = p.network.M,
     bm_sp = round.(bm.species, digits = n_digits),
     bm_total = bm.total,
     bm_pred = sum(bm.species[troph_class.top_predators]),
     bm_cons = sum(bm.species[troph_class.intermediate_consumers]),
     bm_prod = sum(bm.species[troph_class.producers]),
     species_alive = species[troph.alive_species],
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
        species_extirpated = missing,
        tlvl_extirpated = missing,
        species_introduced = missing,
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
           species_extirpated = species_extirpated,
           species_introduced = species_introduced,
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

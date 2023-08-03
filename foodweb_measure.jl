"""
    empirical_interaction_strength(solution, params; kwargs...)

Computes the empirical trophic interaction strength over the `last` timesteps.
It returns the mean, max, min, standard deviation of each trophic interaction in the
network, and return the timeseries as well.

"""
function empirical_interaction_strength(solution, params::ModelParameters; kwargs...)

    measure_on = extract_last_timesteps(solution; kwargs...)

    S = richness(params.network)
    ntimestep = size(measure_on, 2)
    out = zeros(S, S, ntimestep)
    for i in 1:ntimestep
        out[:, :, i] = empirical_interaction_strength(measure_on[:, i], params)
    end

    (
     mean = mean(out, dims = 3)[:, :, 1],
     max = maximum(out, dims = 3)[:, :, 1],
     min = minimum(out, dims = 3)[:, :, 1],
     std = std(out, dims = 3)[:,:, 1],
     all = out
   )

end

function empirical_interaction_strength(B::Vector{Float64}, params::ModelParameters)

    S = size(params.network.species, 1)
    int = zeros(S, S)
    h = params.functional_response.h
    B0 = params.functional_response.B0
    ω = params.functional_response.ω
    x = params.biorates.x
    y = params.biorates.y
    c = params.functional_response.c

    for i in 1:S
        int[i, :] = [
                     ( x[i] * y[i] * B[i] * ω[i, j] * (B[j])^h ) /
                     ( (B0[i])^h + c[i] * B[i] * (B0[i])^h  + sum(ω[i,:] .* (B .^h)))
                     for j in 1:S
                    ]
    end
    sparse(int)
end

"""
```jldoctest
fw = FoodWeb(nichemodel, 10, C = .3)
# My function non weighted is equivalent to the one of EcologicalNetworks pkg
omnivory(fw.A, weighted = false)
EcologicalNetworks.omnivory(UnipartiteNetwork(fw.A))
# Weighted = true correspond to the original definition of Pauly (1987)
omnivory(fw; weighted = false)
omnivory(fw; weighted = true)

# ModelParameters
p = ModelParameters(fw)
omnivory(p.functional_response.ω)
omnivory(p; weighted = false)
```
"""
function omnivory(A; weighted = true)
    # Convert to Bool if preference matrix
    tlvl = trophic_levels(A .!= 0)
    omnivory = Float64[]
    for i in 1:size(A, 1)
        link_indexes = findall(!=(0), A[i,:])
        prey_tlvl = tlvl[link_indexes]
        # Relative preference:
        if sum(A[i,:]) == 0 # if no prey
            rel_pref = []
            push!(omnivory, 0.0)
            out = 0
        else
            # To vector if sparse vector
            prey_tlvl_var = (prey_tlvl .- mean(prey_tlvl)).^2
            if weighted
                rel_pref = Vector(A[i,link_indexes] ./ sum(A[i,:]))
                out = sum(prey_tlvl_var .* rel_pref)
            else
                out = sum(prey_tlvl_var)
            end
            push!(omnivory, out)
        end
    end
    omnivory
end


function omnivory(fw::FoodWeb; kwargs...)
    (species = fw.species, omnivory = omnivory(fw.A; kwargs...))
end

function omnivory(p::ModelParameters; kwargs...)
    omnivory(p.functional_response.ω; kwargs...)
end

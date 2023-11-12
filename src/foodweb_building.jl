
function try_foodweb(S; C = .1, tol_C = .05, n = 5, kwargs...)
    i = 1
    while i <= n
        println("i = $i")
        fw = try
            FoodWeb(nichemodel, S;
                    C = C, tol_C = tol_C,
                    kwargs...
                   )
        catch
            missing
        end
        if !ismissing(fw)
            break
        end
        i += 1
    end
    return fw
end

function add_remove_link_species(foodweb, species::Int; nb_link = 1)
    A = Matrix(foodweb.A)

    L_sp = Vector(A[species,:])
    # Nb of interaction
    nbL_sp = sum(L_sp)
    # Nb of link to add or remove:
    nbL2change = nb_link - nbL_sp
    abs_nbL2change = abs(nbL2change)

    # Where are interactions?
    idxL = findall(>(0), L_sp)
    # Where are absence of interactions ?
    idxNL = findall(==(0), L_sp)

    newA = deepcopy(A) #deepcopy to avoid changing raw foodweb

    if nbL2change < 0
        # Choose randomly link to remove
        L_remove = sample(idxL, abs_nbL2change)
        # Removed the link
        newA[species, L_remove] .= 0
    else nbL2change > 0
        # Choose randomly link to add
        L2add = sample(idxNL, abs_nbL2change)
        newA[species, L2add] .= 1
    end
    newA
end

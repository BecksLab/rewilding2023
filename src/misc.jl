reshape_array(vec) = reshape(vec, (
                                   round(Int, sqrt(length(vec))),
                                   round(Int, sqrt(length(vec)))
                                  )
                            )

function stoch_gc(gc_thre)
    if rand(Distributions.Uniform(0, 1)) < gc_thre
        println("")
        GC.gc()
        ccall(:malloc_trim, Cvoid, (Cint,), 0)
        GC.safepoint()
    end
end

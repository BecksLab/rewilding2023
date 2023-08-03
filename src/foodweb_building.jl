
function try_foodweb(S; C = .1, tol_C = .05, n = 5, kwargs...)
    fw = missing
    local i = 1

    while all([ismissing(fw), i <= n])
        println("i = $i")
        fw = try FoodWeb(nichemodel,
                         S;
                         C = C,
                         tol_C = tol_C,
                         kwargs...
                        )
        catch
            missing
        end
        i = i + 1
    end
    fw
end

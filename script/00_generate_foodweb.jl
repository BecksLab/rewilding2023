using CSV, StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
include("src/foodweb_building.jl")

nrep = 100
S = [20, 30, 40, 60]
C = [0.05, 0.1, 0.15, 0.2]

########################
#  Generate Food-webs  #
########################

rep = 1:nrep
names = (:rep, :richness, :connectance)
param = map(p -> (;Dict(k => v for (k, v) in zip(names, p))...),
            Iterators.product(rep, S, C)
           )[:]

foodweb = map(p -> (rep = p.rep, S = p.richness, C = p.connectance,
                    fw = try_foodweb(p.richness; C = p.connectance,
                                tol_C = .05,
                                check_cycle = false,
                                check_disconnected = true,
                                n = 5
                               ).A
                  ),
            param
            )
df_fw = DataFrame(foodweb)
# Create a foodweb_id
df_fw[!, :fw_id] = 1:nrow(df_fw)

Arrow.write("data/fw_C_S.arrow", df_fw)

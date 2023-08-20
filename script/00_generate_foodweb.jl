println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

#Flag enables all the workers to start on the project of the current dir
dir = "~/rewilding2023/"
flag = "--project=" * dir
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(ncpu - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -2) cores")

@everywhere import Pkg
@everywhere using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
@everywhere include("src/foodweb_building.jl")

import Random.seed!

seed!(22)

nrep = 100
S = [20, 30, 40, 60]
C = [0.05, 0.1, 0.15, 0.2]

########################
#  Generate Food-webs  #
########################
rep = 1:nrep
name = (:rep, :richness, :connectance)
param = map(p -> (;Dict(k => v for (k, v) in zip(name, p))...),
            Iterators.product(rep, S, C)
           )[:]

foodweb = pmap(p -> (rep = p.rep, S = p.richness, C = p.connectance,
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

Arrow.write(dir * "data/fw_C_S.arrow", df_fw)


############################################################################################
#                                    Other experiments                                     #
############################################################################################

# We stick with:
C = 0.15
S = 60


##################
#  Z experiment  #
##################

# Z experiment
Z = [10, 20, 50, 100, 200, 500, 1000]
selected_fw = filter([:S, :C] => (x, y) -> x == S && y == C, df_fw)

rep = 1:nrep
name = (:rep, :Z)
param = map(p -> (;Dict(k => v for (k, v) in zip(name, p))...),
            Iterators.product(rep, Z)
           )[:]

dfZ = DataFrame(param)
dfZ.fw = [selected_fw.fw[i] for i in dfZ.rep]
dfZ.fw_id = 1:nrow(dfZ)
Arrow.write(dir * "data/fw_Z.arrow", dfZ)

####################
#  Propagule size  #
####################

# Propagule size
B_levels = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75]

rep = 1:nrep
name = (:rep, :introduced_size)
param = map(p -> (;Dict(k => v for (k, v) in zip(name, p))...),
            Iterators.product(rep, B_levels)
           )[:]

dfB = DataFrame(param)
dfB.fw = [selected_fw.fw[i] for i in dfB.rep]
dfB.fw_id = 1:nrow(dfB)
Arrow.write(dir * "data/fw_B.arrow", dfB)

################
#  Generality  #
################

G_levels = [1, 3, 5, 7, 9]

rep = 1:nrep
name = (:rep, :nb_link)
param = map(p -> (;Dict(k => v for (k, v) in zip(name, p))...),
            Iterators.product(rep, G_levels)
           )[:]

dfG = DataFrame(param)
dfG.fw = [selected_fw.fw[i] for i in dfG.rep]
dfG.fw_id = 1:nrow(dfG)
Arrow.write(dir * "data/fw_G.arrow", dfG)

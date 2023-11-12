println("Rundir is $(pwd())")

import Pkg
Pkg.instantiate()
using Distributed, Serialization

ncpu = length(Sys.cpu_info())

#Flag enables all the workers to start on the project of the current dir
dir = expanduser("~/rewilding2023/")
flag = "--project=" * dir
#flag = "--project=."
println("Workers run with flag: $(flag)")
addprocs(15 - 1, exeflags=flag)
#addprocs(5, exeflags=flag)
println("Using $(ncpu -2) cores")

@everywhere import Pkg
@everywhere using StatsBase, DataFrames, Arrow, EcologicalNetworksDynamics
@everywhere include("../src/foodweb_building.jl")

import Random.seed!

seed!(22)

nrep = 100
S = [10, 20, 30, 40, 60]
C = [0.05, 0.1, 0.15, 0.2, .25, 0.3]
Z = [10, 20, 50, 100, 200, 500, 1000]
G_levels = [missing, 1, 3, 5, 7, 9]

########################
#  Generate Food-webs  #
########################
rep = 1:nrep
name = (:rep, :richness, :connectance, :PPMR, :nb_link)
param = map(p -> (;Dict(k => v for (k, v) in zip(name, p))...),
            Iterators.product(rep, S, C, Z, G_levels)
           )[:]

# Filter impossible combination of C/S
limitCS = (
           S = S,
           Cmin = round.([(i - 1)/ i^2 + .01 for i in S], digits = 2),
           Cmax = [.33, .24, .22, .17, .12]
          )
# Select good combinations
goodCSparam_v = [
                 findall(x ->
                         (
                          x.connectance >= limitCS.Cmin[i] &&
                          x.connectance <= limitCS.Cmax[i]) &&
                         x.richness == limitCS.S[i],
                         param
                        )
                 for i in 1:length(limitCS.S)
                ]
goodCSparam_idxs = reduce(vcat, goodCSparam_v)
param = param[goodCSparam_idxs]

foodweb = pmap(p -> (
                     rep = p.rep, S = p.richness, C = p.connectance,
                     Z = p.PPMR, nb_link = p.nb_link,
                     fw = try
                         FoodWeb(nichemodel, p.richness; C = p.connectance,
                                 Z = p.PPMR,
                                 tol_C = .03,
                                 check_cycle = true,
                                 check_disconnected = true).A
                     catch
                         missing
                     end
                  ),
               param
            )
df_fw = DataFrame(foodweb)
# Create a foodweb_id
df_fw[!, :sim_id] = 1:nrow(df_fw)
# Remove missing fw
df_fw = df_fw[[!ismissing(df_fw[i, :fw]) for i in 1:nrow(df_fw)], :]

Arrow.write(joinpath(dir, "data/sim_param.arrow"), df_fw)

# Create a toy parameter dataset for testing
param_toy_v = findall(x -> x.rep == 1, param)
param_toy_idxs = reduce(vcat, param_toy_v)
param_toy = param[param_toy_idxs]
foodweb_toy = pmap(p -> (
                     rep = p.rep, S = p.richness, C = p.connectance,
                     Z = p.PPMR, nb_link = p.nb_link,
                     fw = try
                         FoodWeb(nichemodel, p.richness; C = p.connectance,
                                 Z = p.PPMR,
                                 tol_C = .03,
                                 check_cycle = true,
                                 check_disconnected = true).A
                     catch
                         missing
                     end
                  ),
               param_toy
            )
df_fw_toy = DataFrame(foodweb_toy)
# Create a foodweb_id
df_fw_toy[!, :sim_id] = 1:nrow(df_fw_toy)
# Remove missing fw
df_fw_toy = df_fw_toy[[!ismissing(df_fw_toy[i, :fw]) for i in 1:nrow(df_fw_toy)], :]

Arrow.write(joinpath(dir, "data/sim_param_toy.arrow"), df_fw_toy)

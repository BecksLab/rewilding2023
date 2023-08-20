
import Pkg
Pkg.add(["CSV", "StatsBase", "DataFrames", "Arrow"])
Pkg.develop(url = "../src/EcologicalNetworksDynamics.jl-dev")
Pkg.add(["Random", "Plots", "Distributions", "EcologicalNetworksPlots"])
Pkg.add(["DifferentialEquations"])

Pkg.instantiate()

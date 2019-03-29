using Pkg
Pkg.resolve()
Pkg.add("InformationMeasures")
Pkg.add("PyPlot")
Pkg.add("LightGraphs")
Pkg.add("GraphPlot")
Pkg.add("NetworkInference")

# Hopefully this precompiles the libraries
using NetworkInference
using LightGraphs


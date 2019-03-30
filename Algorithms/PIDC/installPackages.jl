using Pkg
Pkg.resolve()

Pkg.add(PackageSpec(name="InformationMeasures", version = "0.3.0"))

Pkg.add(PackageSpec(name="PyPlot", version = "2.8.0"))

Pkg.add(PackageSpec(name="LightGraphs", version = "1.2.0"))
           
Pkg.add(PackageSpec(name = "GraphPlot", version = "0.3.1"))
           
Pkg.add(PackageSpec(name = "NetworkInference", version = "0.1.0"))

# Hopefully this step precompiles the libraries
using NetworkInference
using LightGraphs


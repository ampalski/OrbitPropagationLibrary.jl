module OrbitPropagationLibrary

using OrbitPropagationLibrarySOFA
using DifferentialEquations
using StaticArrays
using DataFrames

include("Utils.jl")
include("TypeDefs.jl")
include("OrbitalElements.jl")
include("Propagation.jl")

const μ = 3.986e5
# for now, replace this with more accurate version (and maybe unitful version) 
# once numeric is up TODO: 

end

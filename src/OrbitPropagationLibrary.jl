module OrbitPropagationLibrary

using OrbitPropagationLibrarySOFA
using DifferentialEquations
using StaticArrays

include("TypeDefs.jl")
include("Propagation.jl")

end

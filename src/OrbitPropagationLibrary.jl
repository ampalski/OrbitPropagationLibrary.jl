module OrbitPropagationLibrary

using OrbitPropagationLibrarySOFA
using DifferentialEquations
using StaticArrays
using DataFrames
using LinearAlgebra

include("Utils.jl")
include("TypeDefs.jl")
include("OrbitalElements.jl")
include("Propagation.jl")

const μ = 3.986004415e5

end

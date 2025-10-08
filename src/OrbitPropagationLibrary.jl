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
const μMOON = 4902.799
const μSUN = 1.32712428e11
const AU = 149597870.7

end

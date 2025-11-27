module OrbitPropagationLibrary

using OrbitPropagationLibrarySOFA
using DifferentialEquations
using StaticArrays
using DataFrames
using LinearAlgebra
using SparseArrays

include("Utils.jl")
include("TypeDefs.jl")
include("OrbitalElements.jl")
include("Propagation.jl")

const μ = 3.986004415e5
const μMOON = 4902.799
const μSUN = 1.32712428e11
const AU = 149597870.7
const REarth = 6378.1363
const RSun = 695700.0
include("EGM2008_200.jl")
export NormGravityModel_S, NormGravityModel_C

end

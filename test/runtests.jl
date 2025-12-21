using OrbitPropagationLibrary
using OrbitPropagationLibrarySOFA
using Aqua
using Test

const tests = [
    "Propagations/Numerical",
    "Extras/aqua",
]
@testset "OrbitPropagationLibrary.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end

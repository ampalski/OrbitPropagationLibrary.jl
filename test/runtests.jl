using OrbitPropagationLibrary
using OrbitPropagationLibrarySOFA
using Test

const tests = [
    "Propagations/Numerical",
]
@testset "OrbitPropagationLibrary.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end

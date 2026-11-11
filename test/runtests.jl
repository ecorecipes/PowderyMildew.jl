using Test
using PowderyMildew

@testset "PowderyMildew" begin
    include("test_hamelin.jl")
    include("test_leaf_model.jl")
    include("test_spatial.jl")
    include("test_interventions.jl")
    include("test_docstrings.jl")
end

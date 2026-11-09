# Every exported name this package owns has a docstring. A name re-exported from another
# package (solve) is that package's to document.
function owned_undocumented(M; allow = Symbol[])
    return filter(n -> Base.which(M, n) === M && n ∉ allow, Base.Docs.undocumented_names(M))
end

@testset "docstrings on every exported name" begin
    @test isempty(owned_undocumented(PowderyMildew))
end

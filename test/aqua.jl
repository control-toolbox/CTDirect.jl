println("Test: Aqua")
using Aqua
using CTDirect

#@testset "Aqua.jl" begin
Aqua.test_all(
    CTDirect;
    ambiguities=false,
    #project_extras = false, seems ignored
    #stale_deps=(ignore=[:SomePackage],),
    deps_compat=(ignore=[:LinearAlgebra, :Unicode],),
    piracies=true,
)
# do not warn about ambiguities in dependencies
Aqua.test_ambiguities(CTDirect)
#end

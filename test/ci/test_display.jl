println("testing: display of DOCP, DOCPCache and ExaBuildCache")

# Load the double integrator problems
if !isdefined(Main, :double_integrator_minenergy)
    include("../problems/double_integrator.jl")
end

# --- fakes and helpers -----------------------------------------------------
# A scheme that is not one of CTDirect's: no description, shown by its name only.
struct DisplayFakeScheme <: CTDirect.Scheme end

const DISPLAY_SCHEMES = (
    :euler,
    :euler_implicit,
    :trapeze,
    :midpoint,
    :gauss_legendre_2_constant_control,
    :gauss_legendre_3_constant_control,
    :gauss_legendre_2,
    :gauss_legendre_3,
)

display_plain(x) = sprint(show, MIME"text/plain"(), x)
display_compact(x) = sprint(show, x)
display_colored(x) = sprint(show, MIME"text/plain"(), x; context=:color => true)
display_strip_ansi(s) = replace(s, r"\e\[[0-9;]*m" => "")

display_cache(ocp, discretizer) = CTSolvers.DOCP.discretize(ocp, discretizer).cache

# Expected outputs, rebuilt from the DOCP fields.
function display_expected_tree(title, docp; scheme_line, steps_line)
    return "$title\n" *
           "├─ scheme: $scheme_line\n" *
           "├─ steps: $steps_line\n" *
           "└─ NLP: $(docp.dim_NLP_variables) variables, $(docp.dim_NLP_constraints) constraints"
end

@testset verbose = true showtiming = true ":display" begin
    ocp = double_integrator_minenergy().ocp

    # --- exact output on the reference case (issue control-toolbox/CTSolvers.jl#337) ---
    @testset "trapeze, grid_size=100" begin
        cache = display_cache(ocp, CTDirect.Collocation(; grid_size=100, scheme=:trapeze))
        @test cache.docp.dim_NLP_variables == 303
        @test cache.docp.dim_NLP_constraints == 204
        tree(title) =
            "$title\n" *
            "├─ scheme: Trapeze (Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable)\n" *
            "├─ steps: 100\n" *
            "└─ NLP: 303 variables, 204 constraints"
        @test display_plain(cache) == tree("DOCPCache")
        @test display_plain(cache.docp) == tree("DOCP")
        @test display_compact(cache) == "DOCPCache(Trapeze, steps=100, nvar=303, ncon=204)"
        @test display_compact(cache.docp) == "DOCP(Trapeze, steps=100, nvar=303, ncon=204)"
        @test sprint(print, cache) == display_compact(cache)
        @test repr(cache.docp) == display_compact(cache.docp)
    end

    # --- every collocation scheme ----------------------------------------------
    @testset "scheme $scheme" for scheme in DISPLAY_SCHEMES
        cache = display_cache(ocp, CTDirect.Collocation(; grid_size=10, scheme=scheme))
        docp = cache.docp
        name = string(nameof(typeof(docp.discretization)))
        description = CTDirect._scheme_description(docp.discretization)
        @test description isa String
        @test description == docp.discretization.info
        @test display_plain(cache) == display_expected_tree(
            "DOCPCache", docp; scheme_line="$name ($description)", steps_line="10"
        )
        @test display_compact(cache) ==
            "DOCPCache($name, steps=10, nvar=$(docp.dim_NLP_variables), " *
              "ncon=$(docp.dim_NLP_constraints))"
    end

    # --- direct shooting: controls per step ------------------------------------
    @testset "direct shooting" begin
        cache = display_cache(ocp, CTDirect.DirectShooting(; grid_size=20, control_steps=3))
        docp = cache.docp
        @test docp.time.control_steps == 3
        @test occursin("├─ steps: 20 (3 controls per step)\n", display_plain(cache))
        @test occursin(", steps=20, control_steps=3, ", display_compact(cache))

        # One control per step (the default) is not mentioned.
        cache = display_cache(ocp, CTDirect.DirectShooting(; grid_size=20))
        @test occursin("├─ steps: 20\n", display_plain(cache))
        @test !occursin("control", display_plain(cache))
        @test !occursin("control_steps", display_compact(cache))
    end

    # --- a scheme without description ------------------------------------------
    @testset "scheme description fallback" begin
        @test CTDirect._scheme_description(DisplayFakeScheme()) === nothing
    end

    # --- ExaBuildCache ---------------------------------------------------------
    @testset "ExaBuildCache" begin
        docp = CTSolvers.DOCP.discretize(ocp, CTDirect.Collocation(; grid_size=50))
        built = CTSolvers.Optimization.build_model(
            docp, CTModels.build_initial_guess(ocp, nothing), CTSolvers.Modelers.Exa()
        )
        @test built.cache isa CTDirect.ExaBuildCache
        @test display_compact(built.cache) == "ExaBuildCache(getter)"
        @test display_plain(built.cache) == "ExaBuildCache(getter)"
    end

    # --- non-regression: size independent of the grid, no type parameters -----
    @testset "size independent of the grid" begin
        small = display_cache(ocp, CTDirect.Collocation(; grid_size=10, scheme=:trapeze))
        large = display_cache(ocp, CTDirect.Collocation(; grid_size=1000, scheme=:trapeze))
        # 10 -> 1000 steps: 2 more digits for steps, nvar and ncon
        @test length(display_plain(large)) - length(display_plain(small)) == 6
        @test length(display_compact(large)) - length(display_compact(small)) == 6
        @test length(display_plain(large)) < 200
        for x in (large, large.docp), s in (display_plain(x), display_compact(x))
            @test !occursin('{', s)
            @test !endswith(s, '\n')
        end
    end

    # --- colors ----------------------------------------------------------------
    @testset "colors" begin
        cache = display_cache(ocp, CTDirect.Collocation(; grid_size=10))
        for x in (cache, cache.docp)
            @test !occursin('\e', display_plain(x))
            @test occursin('\e', display_colored(x))
            @test display_strip_ansi(display_colored(x)) == display_plain(x)
            @test display_strip_ansi(sprint(show, x; context=:color => true)) ==
                display_compact(x)
        end
    end

    # --- composition with the CTSolvers display (CTSolvers.jl#341) --------------
    @testset "DiscretizedModel display" begin
        docp = CTSolvers.DOCP.discretize(ocp, CTDirect.Collocation(; grid_size=1000))
        # The short display of DiscretizedModel comes with CTSolvers.jl#341; before it,
        # Julia's default display prints the whole cache.
        has_display =
            which(show, (IO, MIME"text/plain", typeof(docp))).module === CTSolvers.DOCP
        if has_display
            @test length(display_plain(docp)) < 300
        else
            @test_skip length(display_plain(docp)) < 300
        end
    end
end

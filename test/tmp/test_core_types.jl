# Unit tests for CTDirect core types (integrator schemes and Collocation discretizer).
function test_ctdirect_core_types()

    # ========================================================================
    # TYPE HIERARCHY
    # ========================================================================

    Test.@testset "type hierarchy" verbose=VERBOSE showtiming=SHOWTIMING begin
        # AbstractIntegratorScheme should be abstract
        Test.@test isabstracttype(CTDirect.AbstractIntegratorScheme)

        # Concrete schemes should be subtypes
        Test.@test CTDirect.MidpointScheme <: CTDirect.AbstractIntegratorScheme
        Test.@test CTDirect.TrapezoidalScheme <: CTDirect.AbstractIntegratorScheme

        # Trapeze is an alias to Trapezoidal
        Test.@test CTDirect.TrapezeScheme === CTDirect.TrapezoidalScheme

        # AbstractOptimalControlDiscretizer should be abstract
        Test.@test isabstracttype(CTDirect.AbstractOptimalControlDiscretizer)

        # Collocation should be a concrete discretizer subtype
        Test.@test CTDirect.Collocation <: CTDirect.AbstractOptimalControlDiscretizer
    end

    # ========================================================================
    # COLLOCATION BEHAVIOUR
    # ========================================================================

    Test.@testset "Collocation options and scheme_symbol" verbose=VERBOSE showtiming=SHOWTIMING begin
        # Build a Collocation and read its default options via the generic
        # options API. This keeps the test aligned with the public access
        # pattern instead of calling low-level helpers directly.
        default_colloc = CTDirect.Collocation()
        default_grid = CTModels.get_option_value(default_colloc, :grid)
        default_scheme = CTModels.get_option_value(default_colloc, :scheme)
        default_l2m = CTModels.get_option_value(default_colloc, :lagrange_to_mayer)

        # Sanity checks on defaults
        Test.@test default_grid isa Int
        Test.@test default_grid > 0
        Test.@test default_scheme isa CTDirect.AbstractIntegratorScheme
        Test.@test default_scheme isa CTDirect.MidpointScheme
        Test.@test default_l2m === false

        # Explicitly construct Collocation with given grid and scheme
        colloc = CTDirect.Collocation(;
            grid=default_grid, scheme=default_scheme, lagrange_to_mayer=true
        )

        # Collocation options should expose the stored grid and scheme via options_values
        Test.@test CTModels.get_option_value(colloc, :grid) == default_grid
        Test.@test CTModels.get_option_value(colloc, :scheme) === default_scheme
        Test.@test CTModels.get_option_value(colloc, :lagrange_to_mayer) === true

        # The type parameter of Collocation should reflect the concrete scheme type
        Test.@test default_colloc isa CTDirect.Collocation{CTDirect.MidpointScheme}
        Test.@test colloc isa CTDirect.Collocation{CTDirect.MidpointScheme}
    end

    Test.@testset "discretizer symbols and registry" verbose=VERBOSE showtiming=SHOWTIMING begin
        # get_symbol should return :collocation for the Collocation type and instances.
        Test.@test CTModels.get_symbol(CTDirect.Collocation) == :collocation
        Test.@test CTModels.get_symbol(CTDirect.Collocation()) == :collocation

        # The registered discretizer types should include Collocation.
        regs = CTDirect.registered_discretizer_types()
        Test.@test CTDirect.Collocation in regs

        syms = CTDirect.discretizer_symbols()
        Test.@test :collocation in syms

        # build_discretizer_from_symbol should construct a Collocation
        # discretizer. Use the defaults read from a Collocation instance so
        # that we stay on the generic options API.
        base_disc = CTDirect.Collocation()
        default_grid = CTModels.get_option_value(base_disc, :grid)
        default_scheme = CTModels.get_option_value(base_disc, :scheme)
        disc = CTDirect.build_discretizer_from_symbol(
            :collocation; grid=default_grid, scheme=default_scheme
        )
        Test.@test disc isa CTDirect.Collocation
        Test.@test CTModels.get_option_value(disc, :grid) == default_grid
        Test.@test CTModels.get_option_value(disc, :scheme) === default_scheme
    end

    Test.@testset "build_discretizer_from_symbol unknown symbol" verbose=VERBOSE showtiming=SHOWTIMING begin
        err = nothing
        try
            CTDirect.build_discretizer_from_symbol(:foo)
        catch e
            err = e
        end
        Test.@test err isa CTBase.IncorrectArgument

        buf = sprint(showerror, err)
        Test.@test occursin("Unknown discretizer symbol", buf)
        Test.@test occursin("foo", buf)
        Test.@test occursin("collocation", buf)
    end

    Test.@testset "Collocation default_options and option_default" verbose=VERBOSE showtiming=SHOWTIMING begin
        opts = CTModels.default_options(CTDirect.Collocation)

        # Read the defaults through the generic options API on a default
        # Collocation instance instead of calling low-level helpers.
        base_disc = CTDirect.Collocation()
        default_grid = CTModels.get_option_value(base_disc, :grid)
        default_scheme = CTModels.get_option_value(base_disc, :scheme)
        default_l2m = CTModels.get_option_value(base_disc, :lagrange_to_mayer)

        Test.@test opts.grid == default_grid
        Test.@test opts.scheme === default_scheme
        Test.@test opts.lagrange_to_mayer === default_l2m

        # Type-based and instance-based views of the options metadata should agree.
        colloc_type = typeof(base_disc)

        opts_from_type = CTModels.default_options(CTDirect.Collocation)
        opts_from_inst = CTModels.default_options(colloc_type)
        Test.@test opts_from_inst == opts_from_type

        keys_from_type = CTModels.options_keys(CTDirect.Collocation)
        keys_from_inst = CTModels.options_keys(colloc_type)
        Test.@test Set(keys_from_inst) == Set(keys_from_type)

        Test.@test CTModels.option_default(:grid, CTDirect.Collocation) == default_grid
        Test.@test CTModels.option_default(:scheme, CTDirect.Collocation) ===
            default_scheme
        Test.@test CTModels.option_default(:grid, colloc_type) == default_grid
        Test.@test CTModels.option_default(:scheme, colloc_type) === default_scheme

        Test.@test CTModels.option_default(:lagrange_to_mayer, CTDirect.Collocation) ===
            false
        Test.@test CTModels.option_default(:lagrange_to_mayer, colloc_type) === false
    end
end
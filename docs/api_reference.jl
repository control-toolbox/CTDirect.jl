"""
    generate_api_reference(src_dir::String, ext_dir::String)

Generate the API reference documentation for CTDirect.
Returns the list of pages.

Keep the file lists below in sync with `src/` when files are added, removed, or renamed:
only symbols whose source file is listed get documented (a missing path is silently
ignored, not an error).
"""
function generate_api_reference(src_dir::String, ext_dir::String)
    # Helper to build absolute paths
    src(files...) = [abspath(joinpath(src_dir, f)) for f in files]
    ext(files...) = [abspath(joinpath(ext_dir, f)) for f in files]

    # Symbols to exclude from documentation
    EXCLUDE_SYMBOLS = Symbol[
        :include,
        :eval,
    ]

    pages = [

        # ───────────────────────────────────────────────────────────────────
        # Core Module
        # ───────────────────────────────────────────────────────────────────
        CTBase.automatic_reference_documentation(;
            subdirectory="api",
            primary_modules=[
                CTDirect => src(
                    "CTDirect.jl",
                    "DOCP_cache.jl",
                    "DOCP_data.jl",
                    "DOCP_functions.jl",
                    "DOCP_variables.jl",
                ),
            ],
            exclude=EXCLUDE_SYMBOLS,
            public=false,
            private=true,
            title="Core Module",
            title_in_menu="Core",
            filename="core",
        ),

        # ───────────────────────────────────────────────────────────────────
        # Collocation
        # ───────────────────────────────────────────────────────────────────
        CTBase.automatic_reference_documentation(;
            subdirectory="api",
            primary_modules=[
                CTDirect => src(
                    "collocation.jl",
                    "direct_shooting.jl",
                ),
            ],
            exclude=EXCLUDE_SYMBOLS,
            public=false,
            private=true,
            title="Collocation Methods",
            title_in_menu="Collocation",
            filename="collocation",
        ),

        # ───────────────────────────────────────────────────────────────────
        # Discretization Schemes
        # ───────────────────────────────────────────────────────────────────
        CTBase.automatic_reference_documentation(;
            subdirectory="api",
            primary_modules=[
                CTDirect => src(
                    joinpath("ode", "common.jl"),
                    joinpath("ode", "euler.jl"),
                    joinpath("ode", "irk.jl"),
                    joinpath("ode", "irk_stagewise.jl"),
                    joinpath("ode", "midpoint.jl"),
                    joinpath("ode", "trapeze.jl"),
                ),
            ],
            exclude=EXCLUDE_SYMBOLS,
            public=false,
            private=true,
            title="Discretization Schemes",
            title_in_menu="Discretization",
            filename="discretization",
        ),

    ]

    return pages
end

"""
    with_api_reference(f::Function, src_dir::String, ext_dir::String)

Generates the API reference, executes `f(pages)`, and cleans up generated files.
"""
function with_api_reference(f::Function, src_dir::String, ext_dir::String)
    pages = generate_api_reference(src_dir, ext_dir)
    try
        f(pages)
    finally
        # Clean up generated files
        docs_src = abspath(joinpath(@__DIR__, "src"))
        _cleanup_pages(docs_src, pages)
    end
end

function _cleanup_pages(docs_src::String, pages)
    for p in pages
        val = last(p)
        if val isa AbstractString
            fname = endswith(val, ".md") ? val : val * ".md"
            full_path = joinpath(docs_src, fname)
            if isfile(full_path)
                rm(full_path)
                println("Removed temporary API doc: $full_path")
            end
        elseif val isa AbstractVector
            _cleanup_pages(docs_src, val)
        end
    end
end

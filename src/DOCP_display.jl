# DOCP display
#
# `Base.show` methods for `DOCP`, `DOCPCache` and `ExaBuildCache`. The output never
# contains type parameters, bounds or time grids, so its length does not depend on the
# grid size.

"""
$(TYPEDSIGNATURES)

Description of a discretization scheme (order, properties), or `nothing` for a scheme
that does not provide one.

The schemes of CTDirect store it in their `info` field; any other subtype of
[`CTDirect.Scheme`](@ref) falls back to `nothing` and is displayed by its type name only.

# Arguments
- `scheme::Scheme`: The discretization scheme.

# Returns
- `Union{String,Nothing}`: The description, or `nothing`.

# Example
```julia-repl
julia> CTDirect._scheme_description(docp.discretization)
"Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable"
```
"""
_scheme_description(::Scheme) = nothing
function _scheme_description(
    scheme::Union{Euler,Midpoint,Trapeze,GenericIRK,GenericIRKStagewise}
)
    return scheme.info
end

"""
$(TYPEDSIGNATURES)

Print the number of time steps of `docp`, followed by the number of controls per step
when it is greater than one (direct shooting).

# Arguments
- `io::IO`: Output stream.
- `docp::DOCP`: The discretized problem.

# Returns
- `Nothing`: Prints to `io`.
"""
function _print_steps(io::IO, docp::DOCP)
    print(io, docp.time.steps)
    docp.time.control_steps > 1 &&
        print(io, " (", docp.time.control_steps, " controls per step)")
    return nothing
end

"""
$(TYPEDSIGNATURES)

Print `title` followed by the tree of `docp`: the scheme with its description, the
number of steps, and the numbers of NLP variables and constraints.

Colors follow the active CTBase palette and are emitted only when `io` has `:color`
set. No trailing newline is printed.

# Arguments
- `io::IO`: Output stream.
- `title::AbstractString`: First line of the tree.
- `docp::DOCP`: The discretized problem.

# Returns
- `Nothing`: Prints to `io`.

# Example
```julia-repl
julia> CTDirect._show_docp_tree(stdout, "DOCP", docp)
DOCP
├─ scheme: Trapeze (Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable)
├─ steps: 100
└─ NLP: 303 variables, 204 constraints
```
"""
function _show_docp_tree(io::IO, title::AbstractString, docp::DOCP)
    fmt = Core.get_format_codes(io)
    function branch(prefix, label)
        return print(
            io, '\n', fmt.muted, prefix, fmt.reset, fmt.name, label, fmt.reset, ": "
        )
    end
    print(io, fmt.name, title, fmt.reset)
    branch("├─ ", "scheme")
    print(io, nameof(typeof(docp.discretization)))
    description = _scheme_description(docp.discretization)
    description === nothing || print(io, " (", description, ")")
    branch("├─ ", "steps")
    _print_steps(io, docp)
    branch("└─ ", "NLP")
    print(
        io, docp.dim_NLP_variables, " variables, ", docp.dim_NLP_constraints, " constraints"
    )
    return nothing
end

"""
$(TYPEDSIGNATURES)

Print `name` followed by the scheme, the number of steps and the NLP sizes of `docp`,
on one line.

# Arguments
- `io::IO`: Output stream.
- `name::AbstractString`: Name printed before the parentheses.
- `docp::DOCP`: The discretized problem.

# Returns
- `Nothing`: Prints to `io`.

# Example
```julia-repl
julia> CTDirect._show_docp_compact(stdout, "DOCP", docp)
DOCP(Trapeze, steps=100, nvar=303, ncon=204)
```
"""
function _show_docp_compact(io::IO, name::AbstractString, docp::DOCP)
    fmt = Core.get_format_codes(io)
    print(io, fmt.name, name, fmt.reset, "(", nameof(typeof(docp.discretization)))
    print(io, ", steps=", docp.time.steps)
    docp.time.control_steps > 1 && print(io, ", control_steps=", docp.time.control_steps)
    print(io, ", nvar=", docp.dim_NLP_variables, ", ncon=", docp.dim_NLP_constraints, ")")
    return nothing
end

"""
$(TYPEDSIGNATURES)

Compact display of a [`CTDirect.DOCP`](@ref): scheme, number of steps and NLP sizes on one line,
without type parameters.

# Example
```julia-repl
julia> print(docp)
DOCP(Trapeze, steps=100, nvar=303, ncon=204)
```
"""
Base.show(io::IO, docp::DOCP) = _show_docp_compact(io, "DOCP", docp)

"""
$(TYPEDSIGNATURES)

Display a [`CTDirect.DOCP`](@ref) as a tree: the scheme with its description, the number of
steps (and of controls per step for direct shooting), and the numbers of NLP variables
and constraints. Bounds and time grids are not shown, so the output does not depend on
the grid size.

# Example
```julia-repl
julia> docp
DOCP
├─ scheme: Trapeze (Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable)
├─ steps: 100
└─ NLP: 303 variables, 204 constraints
```
"""
Base.show(io::IO, ::MIME"text/plain", docp::DOCP) = _show_docp_tree(io, "DOCP", docp)

"""
$(TYPEDSIGNATURES)

Compact display of a [`CTDirect.DOCPCache`](@ref): the summary of its [`CTDirect.DOCP`](@ref) on one line.

# Example
```julia-repl
julia> print(docp.cache)
DOCPCache(Trapeze, steps=100, nvar=303, ncon=204)
```
"""
Base.show(io::IO, cache::DOCPCache) = _show_docp_compact(io, "DOCPCache", cache.docp)

"""
$(TYPEDSIGNATURES)

Display a [`CTDirect.DOCPCache`](@ref) as the tree of its [`CTDirect.DOCP`](@ref).

# Example
```julia-repl
julia> docp.cache
DOCPCache
├─ scheme: Trapeze (Implicit Trapeze aka Crank-Nicolson, 2nd order, A-stable)
├─ steps: 100
└─ NLP: 303 variables, 204 constraints
```
"""
function Base.show(io::IO, ::MIME"text/plain", cache::DOCPCache)
    return _show_docp_tree(io, "DOCPCache", cache.docp)
end

"""
$(TYPEDSIGNATURES)

Display an [`CTDirect.ExaBuildCache`](@ref) on one line. The getter is not shown: its type
carries the whole ExaModels expression tree.

# Example
```julia-repl
julia> built.cache
ExaBuildCache(getter)
```
"""
function Base.show(io::IO, ::ExaBuildCache)
    fmt = Core.get_format_codes(io)
    print(io, fmt.name, "ExaBuildCache", fmt.reset, "(getter)")
    return nothing
end

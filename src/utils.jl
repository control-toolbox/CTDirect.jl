"""
$(TYPEDSIGNATURES)

Return the version number of the current package as a string.

# Returns

- `version::String`: The version of the current module.

# Example

```julia-repl
julia> version()
"1.2.3"
```
"""
version() = string(pkgversion(@__MODULE__))

"""
$(TYPEDSIGNATURES)

Check whether a package with the given name is currently loaded.

# Arguments

- `pkg::String`: The name of the package to check.

# Returns

- `is_loaded::Bool`: `true` if the package is loaded, `false` otherwise.

# Example

```julia-repl
julia> package_loaded("LinearAlgebra")
true
```
"""
package_loaded(pkg::String) = any(k -> k.name == pkg, keys(Base.loaded_modules))

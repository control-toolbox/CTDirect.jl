version() = string(pkgversion(@__MODULE__))

package_loaded(pkg::String) = any(k -> k.name == pkg, keys(Base.loaded_modules))
# CTDirect API

## Index

```@index
Pages   = ["api.md"]
Modules = [CTDirect]
Order = [:module, :constant, :type, :function, :macro]
```

## Available methods

```@example
using CTDirect
available_methods()
```

## Documentation

```@autodocs
Modules = [CTDirect]
Order = [:module, :constant, :type, :function, :macro]
Private = false
```

<! -- manually add docstrings from package extensions -->
<! -- does not work for solve which is not in CTDirect -->
```@docs
solve
save_OCP_solution
load_OCP_solution
export_OCP_solution
read_OCP_solution
```

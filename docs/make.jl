using Documenter
using CTDirect
using CTBase

DocMeta.setdocmeta!(CTBase, :DocTestSetup, :(using CTBase); recursive = true)
DocMeta.setdocmeta!(CTDirect, :DocTestSetup, :(using CTDirect); recursive = true)

makedocs(
    warnonly = [:cross_references, :autodocs_block],
    sitename = "CTDirect.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        #"Introduction"  => "index.md",
        #"Methods and Options" => ["rk.md",
        #                          "optimization.md",
        #                          "solve-options.md"],
        #"API"           => "api.md",
        "Tutorial"      => "tutorial.md",
        #"Developers"    => "dev-api.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/CTDirect.jl.git",
    devbranch = "main"
)

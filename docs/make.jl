using Documenter
using CTDirect
using CTBase

DocMeta.setdocmeta!(CTBase, :DocTestSetup, :(using CTBase); recursive = true)
DocMeta.setdocmeta!(CTDirect, :DocTestSetup, :(using CTDirect); recursive = true)

makedocs(
    sitename = "CTDirect.jl",
    format = Documenter.HTML(prettyurls = false),
    size_threshold = nothing,
    pages = [
        "Introduction"  => "index.md",
        #"Methods and Options" => ["rk.md",
        #                          "optimization.md",
        #                          "solve-options.md"],
        "API"           => "api.md",
        "Tutorial"      => "tutorial.md",
        #"Developers"    => "dev-api.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/CTDirect.jl.git",
    devbranch = "main"
)

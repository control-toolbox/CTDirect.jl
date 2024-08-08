using Documenter
using DocumenterMermaid
using CTDirect
using CTBase

using NLPModelsIpopt
using HSL
using JLD2
using JSON3
using Plots

# to add docstrings from external packages
DocMeta.setdocmeta!(CTBase, :DocTestSetup, :(using CTBase); recursive = true)
DocMeta.setdocmeta!(CTDirect, :DocTestSetup, :(using CTDirect); recursive = true)

makedocs(
    warnonly = [:cross_references, :autodocs_block],
    sitename = "CTDirect.jl",
    format = Documenter.HTML(
        prettyurls = false,
        size_threshold_ignore = ["api-ctbase.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages = [
        "Introduction"  => "index.md",
        "API"           => "api.md",
        "Developers"    => "dev-api.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/CTDirect.jl.git",
    devbranch = "main"
)

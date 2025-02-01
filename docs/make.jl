using Documenter
using DocumenterMermaid
using CTDirect
using CTBase

using NLPModelsIpopt
using HSL

# to add docstrings from external packages
DocMeta.setdocmeta!(CTBase, :DocTestSetup, :(using CTBase); recursive=true)
DocMeta.setdocmeta!(CTDirect, :DocTestSetup, :(using CTDirect); recursive=true)

repo_url = "github.com/control-toolbox/CTDirect.jl"

makedocs(;
    warnonly=[:cross_references, :autodocs_block],
    sitename="CTDirect.jl",
    format=Documenter.HTML(
        repolink="https://" * repo_url,
        prettyurls=false,
        size_threshold_ignore=["api-ctbase.md"],
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=["Introduction" => "index.md", "API" => "api.md", "Developers" => "dev-api.md"],
)

deploydocs(;
    repo=repo_url * ".git", devbranch="main"
)

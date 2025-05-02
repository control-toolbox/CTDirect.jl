using Documenter
using DocumenterMermaid

using CTDirect

# to add docstrings from external packages
DocMeta.setdocmeta!(CTDirect, :DocTestSetup, :(using CTDirect); recursive=true)

repo_url = "github.com/control-toolbox/CTDirect.jl"

makedocs(;
    warnonly=[:cross_references, :autodocs_block],
    sitename="CTDirect.jl",
    format=Documenter.HTML(
        repolink="https://"*repo_url,
        prettyurls=false,
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=["Introduction" => "index.md", "API" => "api.md", "Developers" => "dev-api.md"],
)

deploydocs(;
    repo=repo_url*".git", devbranch="main"
)

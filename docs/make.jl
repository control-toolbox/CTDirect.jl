using Documenter
using DocumenterMermaid
using CTDirect

repo_url = "github.com/control-toolbox/CTDirect.jl"

API_PAGES = [
    "core_types.md",      
    "collocation_core.md",    
    "collocation.md",
    "common.md",
    "docp.md",
    "euler.md",
    "irk.md",
    "midpoint.md",
    "trapeze.md",
]

makedocs(;
    warnonly=[:cross_references, :autodocs_block],
    sitename="CTDirect.jl",
    format=Documenter.HTML(;
        repolink="https://" * repo_url,
        prettyurls=false,
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=["Introduction" => "index.md", "API" => API_PAGES],
)

deploydocs(; repo=repo_url * ".git", devbranch="main")

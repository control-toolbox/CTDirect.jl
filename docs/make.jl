using Documenter
using DocumenterMermaid
using CTDirect
using NLPModelsIpopt
using NLPModelsKnitro
using MadNLP

# to add docstrings from external packages
const CTSolveExtIpopt = Base.get_extension(CTDirect, :CTSolveExtIpopt)
const CTSolveExtKnitro = Base.get_extension(CTDirect, :CTSolveExtKnitro)
const CTSolveExtMadNLP = Base.get_extension(CTDirect, :CTSolveExtMadNLP)
Modules = [
    CTSolveExtIpopt,
    CTSolveExtKnitro,
    CTSolveExtMadNLP,
]
for Module in Modules
    isnothing(DocMeta.getdocmeta(Module, :DocTestSetup)) &&
        DocMeta.setdocmeta!(Module, :DocTestSetup, :(using $Module); recursive=true)
end

repo_url = "github.com/control-toolbox/CTDirect.jl"

API_PAGES = [
    "common.md",
    "ctsolveext_ipopt.md",
    "ctsolveext_knitro.md",
    "ctsolveext_madnlp.md",
    "default.md",
    "docp.md",
    "euler.md",
    "irk.md",
    "midpoint.md",
    "solution.md",
    "solve.md",
    "trapeze.md",
]

makedocs(;
    warnonly=[:cross_references, :autodocs_block],
    sitename="CTDirect.jl",
    format=Documenter.HTML(
        repolink="https://" * repo_url,
        prettyurls=false,
        assets=[
            asset("https://control-toolbox.org/assets/css/documentation.css"),
            asset("https://control-toolbox.org/assets/js/documentation.js"),
        ],
    ),
    pages=["Introduction" => "index.md", "API" => API_PAGES],
)

deploydocs(;
    repo=repo_url * ".git", devbranch="main"
)

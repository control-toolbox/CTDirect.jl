using Documenter
using DocumenterMermaid
using CTDirect
using ADNLPModels
using ExaModels
using NLPModelsIpopt
using NLPModelsKnitro
using MadNLP

# to add docstrings from external packages
const CTDirectExtADNLP = Base.get_extension(CTDirect, :CTDirectExtADNLP)
const CTDirectExtExa = Base.get_extension(CTDirect, :CTDirectExtExa)
const CTDirectExtIpopt = Base.get_extension(CTDirect, :CTDirectExtIpopt)
const CTDirectExtKnitro = Base.get_extension(CTDirect, :CTDirectExtKnitro)
const CTDirectExtMadNLP = Base.get_extension(CTDirect, :CTDirectExtMadNLP)
Modules = [
    CTDirectExtADNLP, CTDirectExtExa, CTDirectExtIpopt, CTDirectExtKnitro, CTDirectExtMadNLP
]
for Module in Modules
    isnothing(DocMeta.getdocmeta(Module, :DocTestSetup)) &&
        DocMeta.setdocmeta!(Module, :DocTestSetup, :(using $Module); recursive=true)
end

repo_url = "github.com/control-toolbox/CTDirect.jl"

API_PAGES = [
    "common.md",
    "ctdirectext_adnlp.md",
    "ctdirectext_exa.md",
    "ctdirectext_ipopt.md",
    "ctdirectext_knitro.md",
    "ctdirectext_madnlp.md",
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

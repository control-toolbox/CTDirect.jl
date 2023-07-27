using Documenter
using CTDirect

makedocs(
    sitename = "CTDirect.jl",
    format = Documenter.HTML(prettyurls = false),
    pages = [
        "Introduction"  => "index.md",
        "Methods and Options" => ["rk.md",
                                  "optimization.md",
                                  "solve-options.md"],
        "API"           => "api.md",
        #"Developers"    => "dev-api.md",
    ]
)

deploydocs(
    repo = "github.com/control-toolbox/CTDirect.jl.git",
    devbranch = "main"
)

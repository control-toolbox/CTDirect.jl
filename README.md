# CTDirect.jl

[ci-img]: https://github.com/control-toolbox/CTDirect.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/control-toolbox/CTDirect.jl/actions/workflows/CI.yml?query=branch%3Amain

[co-img]: https://codecov.io/gh/control-toolbox/CTDirect.jl/branch/main/graph/badge.svg?token=6J4SJL2SFG
[co-url]: https://codecov.io/gh/control-toolbox/CTDirect.jl

[doc-dev-img]: https://img.shields.io/badge/docs-dev-8A2BE2.svg
[doc-dev-url]: https://control-toolbox.org/CTDirect.jl/dev/

[doc-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[doc-stable-url]: https://control-toolbox.org/CTDirect.jl/stable/

[release-img]: https://juliahub.com/docs/General/CTDirect/stable/version.svg
[release-url]: https://github.com/control-toolbox/CTDirect.jl/releases

[pkg-eval-img]: https://img.shields.io/badge/Julia-package-purple
[pkg-eval-url]: https://juliahub.com/ui/Packages/General/CTDirect

[deps-img]: https://juliahub.com/docs/General/CTDirect/stable/deps.svg
[deps-url]: https://juliahub.com/ui/Packages/General/CTDirect?t=2

[licence-img]: https://img.shields.io/badge/License-MIT-yellow.svg
[licence-url]: https://github.com/control-toolbox/CTDirect.jl/blob/master/LICENSE

The CTDirect.jl package is part of the [control-toolbox ecosystem](https://github.com/control-toolbox).
The control-toolbox ecosystem gathers Julia packages for mathematical control and applications. The root package is [OptimalControl.jl](https://github.com/control-toolbox/OptimalControl.jl) which aims to provide tools to modelise and solve optimal control problems with ordinary differential equations by direct and indirect methods.

[![doc OptimalControl.jl](https://img.shields.io/badge/Documentation-OptimalControl.jl-blue)](http://control-toolbox.org/OptimalControl.jl)

| **Name**          | **Badge**         |
:-------------------|:------------------|
| Documentation     | [![Documentation][doc-stable-img]][doc-stable-url] [![Documentation][doc-dev-img]][doc-dev-url]                   | 
| Code Status       | [![Build Status][ci-img]][ci-url] [![Covering Status][co-img]][co-url] [![pkgeval][pkg-eval-img]][pkg-eval-url]  |
| Dependencies      | [![deps][deps-img]][deps-url] |
| Licence           | [![License: MIT][licence-img]][licence-url]   |
| Release           | [![Release][release-img]][release-url]        |

## Installation

To install CTDirect.jl please 
<a href="https://docs.julialang.org/en/v1/manual/getting-started/">open Julia's interactive session (known as REPL)</a> 
and press <kbd>]</kbd> key in the REPL to use the package mode, then add the package:

```julia
julia> ]
pkg> add CTDirect
```

## Contributing

[issue-url]: https://github.com/control-toolbox/CTDirect.jl/issues
[first-good-issue-url]: https://github.com/control-toolbox/CTDirect.jl/contribute

If you think you found a bug or if you have a feature request / suggestion, feel free to open an [issue][issue-url].
Before opening a pull request, please start an issue or a discussion on the topic. 

Contributions are welcomed, check out [how to contribute to a Github project](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project). 
If it is your first contribution, you can also check [this first contribution tutorial](https://github.com/firstcontributions/first-contributions).
You can find first good issues (if any 🙂) [here][first-good-issue-url]. You may find other packages to contribute to at the [control-toolbox organization](https://github.com/control-toolbox).

If you want to ask a question, feel free to start a discussion [here](https://github.com/orgs/control-toolbox/discussions). This forum is for general discussion about this repository and the [control-toolbox organization](https://github.com/control-toolbox).

>[!NOTE]
> If you want to add an application or a package to the control-toolbox ecosystem, please follow this [set up tutorial](https://github.com/control-toolbox/CTApp.jl/discussions/9).

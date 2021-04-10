# DataEnvelopmentAnalysis.jl
![DataEnvelopmentAnalysis logo](docs/src/assets/wordmark.svg "DataEnvelopmentAnalysis logo")

A Julia package for efficiency and productivity measurement using Data Envelopment Analysis (DEA)

| Documentation | Build Status      | Coverage    | Zenodo      |
|:-------------:|:-----------------:|:-----------:|:-----------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] |  [![][githubci-img]][githubci-url] | [![][codecov-img]][codecov-url] | [![][zenodo-img]][zenodo-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://javierbarbero.github.io/DataEnvelopmentAnalysis.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://javierbarbero.github.io/DataEnvelopmentAnalysis.jl/dev

[githubci-img]: https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl/workflows/CI/badge.svg?branch=master
[githubci-url]: https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl/actions?query=workflow%3ACI

[codecov-img]: https://codecov.io/gh/javierbarbero/DataEnvelopmentAnalysis.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/javierbarbero/DataEnvelopmentAnalysis.jl

[zenodo-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.4412261.svg
[zenodo-url]: https://doi.org/10.5281/zenodo.4412261

DataEnvelopmentAnalysis.jl is a Julia package that provides functions for efficiency and productivity measurement using Data Envelopment Analysis (DEA). Particularly, it implements a variety of technical efficiency models, economic efficiency models and productivity change models.

The package is being developed for Julia `1.0` and above on Linux, macOS, and Windows.

The packes uses internally the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modelling language for mathematicall optimization with solvers [GLPK](http://www.gnu.org/software/glpk/) and [Ipopt](https://coin-or.github.io/Ipopt/).

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("DataEnvelopmentAnalysis")
```

## Available models

**Technical efficiency DEA models**

* Radial input and output oriented model.
* Directional distance function model.
* Additive models: weighted additive model, measure of inefficiency proportions (MIP), normalized weighted additive model, range adjusted measure (RAM), bounded adjusted measure (BAM).
* Generalized distance function model.
* Russell graph and oriented model.
* Enhanced Russell Graph Slack Based Measure.
* Modified directional distance function.
- Hölder distance function.
- Reverse directional distance function.

**Economic efficiency DEA models**

* Cost model.
* Revenue model.
* Profit model.
* Profitability model.

**Productivity change models**

* Mamlmquist index.

## Authors

DataEnvelopmentAnalysis.jl is being developed by [Javier Barbero](http://www.javierbarbero.net) and [José Luís Zofío](http://www.joselzofio.net).

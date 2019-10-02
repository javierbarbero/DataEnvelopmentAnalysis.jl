# DataEnvelopmentAnalysis.jl
A Julia package for efficiency and productivity measurement using Data Envelopment Analysis (DEA)

| Documentation | Build Status      | Coverage    |
|:-------------:|:-----------------:|:-----------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] |  [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] | [![][coveralls-img]][coveralls-url] [![][codecov-img]][codecov-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://javierbarbero.github.io/DataEnvelopmentAnalysis.jl/stable

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://javierbarbero.github.io/DataEnvelopmentAnalysis.jl/dev

[travis-img]: https://travis-ci.org/javierbarbero/DataEnvelopmentAnalysis.jl.svg?branch=master
[travis-url]: https://travis-ci.org/javierbarbero/DataEnvelopmentAnalysis.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/ql1ynl90hcixdka7?svg=true
[appveyor-url]: https://ci.appveyor.com/project/javierbarbero/dataenvelopmentanalysis-jl

[coveralls-img]: https://coveralls.io/repos/github/javierbarbero/DataEnvelopmentAnalysis.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/javierbarbero/DataEnvelopmentAnalysis.jl?branch=master

[codecov-img]: https://codecov.io/gh/javierbarbero/DataEnvelopmentAnalysis.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/javierbarbero/DataEnvelopmentAnalysis.jl

DataEnvelopmentAnalysis.jl is a Julia package that provides functions for efficiency and productivity measurement using Data Envelopment Analysis (DEA). Particularly, it implements a variety of technical efficiency models, economic efficiency models and productivity change models.

The package is being developed for Julia `1.0` and above on Linux, macOS, and Windows.

The packes uses internally the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modelling language for mathematicall optimization with solvers [GLPK](http://www.gnu.org/software/glpk/) and [Ipopt](https://coin-or.github.io/Ipopt/). 

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("DataEnvelopmentAnalysis.jl")
```

## Available models

**Technical efficiency DEA models**

* Radial input and output oriented model.
* Directional distance function model.
* Additive models: weighted additive model, measure of inefficiency proportions (MIP), normalized weighted additive model, range adjusted measure (RAM), bounded adjusted measure (BAM).
* Generalized distance function model.

**Economic efficiency DEA models**

* Cost model.
* Revenue model.
* Profit model.
* Profitability model.

**Productivity change models**

* Mamlmquist index.

## Authors

DataEnvelopmentAnalysis.jl is being developed by [Javier Barbero](http://www.javierbarbero.net) and [José Luís Zofío](http://www.joselzofio.net).


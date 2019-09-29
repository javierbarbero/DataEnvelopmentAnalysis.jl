# DataEnvelopmentAnalysis.jl
A Julia package for efficiency and productivity measurement using Data Envelopment Analysis (DEA)

| Build Statud |
|:-----------------:|
|  [![][travis-img]][travis-url] |

[travis-img]: https://travis-ci.org/javierbarbero/DataEnvelopmentAnalysis.jl.svg?branch=master
[travis-url]: https://travis-ci.org/javierbarbero/DataEnvelopmentAnalysis.jl

DataEnvelopmentAnalysis.jl is a Julia package that provides functions for efficiency and productivity measurement using Data Envelopment Analysis (DEA). Particularly, it implements a variety of technical efficiency models, economic efficiency models and productivity change models.

The package is being developed for Julia `1.0` and above on Linux, macOS, and Windows.

The packes uses internally the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modelling language for mathematicall optimization with solvers [GLPK](http://www.gnu.org/software/glpk/) and [Ipopt](https://coin-or.github.io/Ipopt/). 

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add(PackageSpec(url = "https://github.com/javierbarbero/DataEnvelopmentAnalysis.jl", rev = "master"))
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


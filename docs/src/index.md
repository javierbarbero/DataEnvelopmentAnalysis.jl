# DataEnvelopmentAnalysis Documentation

DataEnvelopmentAnalysis.jl is a Julia package that provides functions for efficiency and productivity measurement using Data Envelopment Analysis (DEA). Particularly, it implements a variety of technical efficiency models, economic efficiency models and productivity change models.

The package is being developed for Julia `1.0` and above on Linux, macOS, and Windows.

The packes uses internally the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modelling language for mathematicall optimization with solvers [GLPK](http://www.gnu.org/software/glpk/) and [Ipopt](https://coin-or.github.io/Ipopt/). 

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("DataEnvelopmentAnalysis.jl")
```

## Available models

Technical efficiency DEA models:
```@contents
Pages = ["technical/radial.md", "technical/directional.md", "technical/additive.md", "technical/generalizeddf.md"]
Depth = 2
```

Economic efficiency DEA models:
```@contents
Pages = ["economic/cost.md", "economic/revenue.md", "economic/profit.md", "economic/profitability.md"]
Depth = 1
```

Productivity change models:
```@contents
Pages = ["productivity/malmquist.md"]
Depth = 1
```

## Authors

DataEnvelopmentAnalysis.jl is being developed by [Javier Barbero](http://www.javierbarbero.net) and [José Luís Zofío](http://www.joselzofio.net).


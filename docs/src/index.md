# DataEnvelopmentAnalysis.jl

DataEnvelopmentAnalysis.jl is a Julia package that provides functions for efficiency and productivity measurement using Data Envelopment Analysis (DEA). Particularly, it implements a variety of technical efficiency models, economic efficiency models and productivity change models.

The package is being developed for Julia `1.6` and above on Linux, macOS, and Windows.

The packes uses internally the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modelling language for mathematicall optimization with solvers [GLPK](https://github.com/jump-dev/GLPK.jl) and [Ipopt](https://github.com/jump-dev/Ipopt.jl).

## Installation

The package can be installed with the Julia package manager:
```julia
julia> using Pkg; Pkg.add("DataEnvelopmentAnalysis")
```

## Tutorial

For a tutorial on how to use the package, check the documentation on the [Radial Input Oriented Model](@ref).

## Available models

Technical efficiency DEA models:
```@contents
Pages = ["technical/radial.md", "technical/radialbigdata.md", "technical/directional.md", "technical/additive.md", "technical/generalizeddf.md", "technical/russell.md", "technical/enhancedrussell.md", "technical/modifiedddf.md", "technical/holder.md", "technical/reverseddf.md", "technical/environmental.md"]
Depth = 2
```

Economic efficiency DEA models:
```@contents
Pages = ["economic/cost.md", "economic/revenue.md", "economic/profit.md", "economic/profitability.md"]
Depth = 1
```

Productivity change models:
```@contents
Pages = ["productivity/malmquist.md", "productivity/malmquistluenberger.md"]
Depth = 1
```

Statistical Analysis:
```@contents
Pages = ["statistical/radialboot.md", "statistical/rtstest.md"]
Depth = 1
```

## Extensions

The [BenchmarkingEconomicEfficiency.jl](https://github.com/javierbarbero/BenchmarkingEconomicEfficiency.jl) package provides an extensive set of functions for economic efficiency measurement using Data Envelopment Analysis (DEA).

## Authors

DataEnvelopmentAnalysis.jl is being developed by [Javier Barbero](http://www.javierbarbero.net) and [José Luís Zofío](http://www.joselzofio.net).

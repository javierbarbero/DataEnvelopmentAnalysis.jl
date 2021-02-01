```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Common functions for economic models

```@docs
efficiency(::AbstractEconomicDEAModel, ::Symbol)
targets(::AbstractEconomicDEAModel, ::Symbol)
normfactor
ismonetary(::AbstractEconomicDEAModel)
```

```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Common functions for technical models

```@docs
nobs(::AbstractDEAModel)
ninputs(::AbstractDEAModel) 
noutputs(::AbstractDEAModel) 
Base.names(::AbstractDEAModel)
efficiency(::AbstractTechnicalDEAModel)
slacks
targets(::AbstractTechnicalDEAModel, ::Symbol)
peers
peersmatrix
ispeer
```

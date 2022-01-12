```@meta
CurrentModule = DataEnvelopmentAnalysis
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
multipliers
rts
```

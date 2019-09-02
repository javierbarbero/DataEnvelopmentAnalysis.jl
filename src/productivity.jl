# This file contains types and structures for productivity DEA model
"""
    AbstractProductivityDEAlModel
An abstract type representing a productivity DEA model.
"""
abstract type AbstractProductivityDEAModel end

"""
    nbos(model::AbstractProductivityDEAModel)
Return number of observations of a productivity DEA model.
# Examples
```jldoctest

```
"""
nobs(model::AbstractProductivityDEAModel) = model.n

"""
    ninputs(model::AbstractProductivityDEAModel)
Return number of inputs of a productivity DEA model.
## Examples
```jldoctest

```
"""
ninputs(model::AbstractProductivityDEAModel) = model.m

"""
    noutputs(model::AbstractProductivityDEAModel)
Return number of outputs of a productivity DEA model.
# Examples
```jldoctest

```
"""
noutputs(model::AbstractProductivityDEAModel) = model.s

"""
    nperiods(model::AbstractProductivityDEAModel)
Return number of time periods of a productivity DEA model.
# Examples
```jldoctest

```
"""
nperiods(model::AbstractProductivityDEAModel) = model.periods

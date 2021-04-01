# This file contains functions for the DEA optimizers
"""
    AbstractDEAOptimizer
An abstract type representing a DEA optimizer.
"""
abstract type AbstractDEAOptimizerModel end

"""
    DEAOptimizer(optimizer; time_limit, silent)
An data structure storing the configuration of a DEA optimizer.

# Optimizer specification:
- `LP`: linear programming default optimizer, GLPK.
- `NLP`: nonlinear programmin default optimizer, Ipopt.
- Any JuMP supported solver.

# Optional Arguments
- `time_limit=:60`: time limit in seconds.
- `silent=:true`: suppress optimizer output.

# Examples
```jldoctest
julia> myoptimizer = DEAOptimizer(GLPK.Optimizer, time_limit = 10, silent = true);
```
"""
struct DEAOptimizer <: AbstractDEAOptimizerModel
    optimizer
    time_limit::Float64
    silent::Bool

    function DEAOptimizer(optimizer; time_limit = 60, silent = true)
        
        if optimizer == :LP
            optimizer =  GLPK.Optimizer
        elseif optimizer == :NLP
            optimizer = Ipopt.Optimizer
        end

        new(optimizer, time_limit, silent)
    end

end

"""
    newdeamodel(DEAOptimizer)
Generate a new JuMP model for DEA with the specified optimizer.

This function is used internally and for packages that want to extend the functionality of this package.

# Examples
```jldoctest
julia> deamodel = newdeamodel(DEAOptimizer(GLPK.Optimizer));
```
"""
function newdeamodel(deaoptimizer::DEAOptimizer)::Model
    deamodel = Model(deaoptimizer.optimizer)

    set_time_limit_sec(deamodel, deaoptimizer.time_limit)

    if deaoptimizer.silent
        set_silent(deamodel)
    end

    return deamodel
end

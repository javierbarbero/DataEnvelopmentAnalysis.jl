```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Hölder Distance Function Models

*Briec (1998)* defined technical inefficiency using Hölder norms.

## Hölder L1

In this example we compute the Hölder L1 DEA model under varible returns to scale:
```@example holder
using DataEnvelopmentAnalysis

X = [2; 4; 8; 12; 6; 14; 14; 9.412];

Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

deaholder(X, Y, l = 1, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example holder
holderl1 = deaholder(X, Y, l = 1, rts = :VRS);
```

```@example holder
efficiency(holderl1)
```

The input or output that determines the projection to the frontier is returned with:
```@example holder
efficiency(holderl1, :min)
```
with inputs and outputs numbered sequentially.

## Hölder L2

!!! warning "Requieres a solver that supports SOS constraints"
    The Hölder L2 model requieres a solver that supports SOS constraints, such as [Gurobi](https://github.com/jump-dev/Gurobi.jl). 

    Solving the model with Ipopt will return invalid results.

## Hölder LInf

In this example we compute the Hölder LInf DEA model under varible returns to scale:
```@example holder
X = [2; 4; 8; 12; 6; 14; 14; 9.412];

Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

deaholder(X, Y, l = Inf, rts = :VRS)
```

### deaholder Function Documentation

```@docs
deaholder
```

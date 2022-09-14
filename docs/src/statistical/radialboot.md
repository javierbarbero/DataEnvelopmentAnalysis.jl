```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Bootstrap Radial Model

The bootstrap radial DEA model (Simar and Wilson, 1998) can be calculated with the `deaboot` function, indicating the number of bootstrap replications in the `nreps` parameter. The other parameters work the same as in the radial DEA model.

A random number generator can be specified in the `rng` parameter for reproducibility.

```@example radialboot
using DataEnvelopmentAnalysis
using StableRNGs

X = [2, 4, 3, 5, 6]
Y = [1, 2, 3, 4, 5]

ioboot = deaboot(X, Y, orient = :Input, rts = :VRS, nreps = 200, rng = StableRNG(1234567))
```

!!! warning "Number of bootstrap replications"
    The example above uses 200 bootstrap replications for illustrative purposes. In practice, at least 1000 replications are recommended.

Bias-corrected efficiency scores are returned with the `efficiency` function:
```@example radialboot
efficiency(ioboot)
```

The bias, calculated as the difference between the reference efficiency score and the bias-corrected efficiency score, is returned with the `bias` function:
```@example radialboot
bias(ioboot)
```

Confidence intervals at the $95\%$, or any other desired level, are calculated with the `confint` function: 
```@example radialboot
confint(ioboot, level = 0.95)
```

The optimal bandwidth computed for the model is returned with the `bandwidth` function:
```@example radialboot
bandwidth(ioboot)
```

### deaboot Function Documentation

```@docs
deaboot
bias(::BootstrapRadialDEAModel)
confint(::BootstrapRadialDEAModel)
bandwidth(::BootstrapRadialDEAModel)
```



```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Returns to Scale (RTS) Test

The non-parametric Returns to Scale test (Simar and Wilson, 2002) is based on the bootstrap radial DEA model. It tests the hypothesis that the technology exhibits constant returns to scale (CRS) against the alternative that it is variable returns to scale (VRS):
```math
\begin{aligned}
  H_0 & \text{: Technology is CRS} \\
  H_1 & \text{: Tecnhology is VRS}
\end{aligned}
```

The test can be performed with the `deartstest` function, indicating the number of bootstrap replications in the `nreps` parameter. A random number generator can be specified in the `rng` parameter for reproducibility.

```@example radialboot
using DataEnvelopmentAnalysis
using StableRNGs

X = [2, 4, 3, 5, 6]
Y = [1, 2, 3, 4, 5]

ioboot = deartstest(X, Y, orient = :Input, nreps = 200, rng = StableRNG(1234567))
```

!!! warning "Number of bootstrap replications"
    The example above uses 200 bootstrap replications for illustrative purposes. In practice, at least 1000 replications are recommended.

We reject the null hypothesis if the estimated scale efficiency is less than the critical value. Alternatively, we can guide our decision using the calculated p-value. In the example above, we do not reject the null hypothesis of constant returns to scale.

### deaboot Function Documentation

```@docs
deartstest
criticalvalue(::DEAReturnsToScaleTest)
```



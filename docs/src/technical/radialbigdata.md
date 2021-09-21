```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Radial Big Data Models

When the number of decision-making units is large, traditional DEA models are slow to solve. Khezrimotlagh, Zhu, Cook, and Toloo (2019), propose a framework that reduces the computational time by finding the set of best practices DMUs from a subsample and evaluating the rest of the decision-making units with respect to the best performers.

The proposed framework includes five steps:
1. Select a subsample of DMU.
2. Find the best practices in the subsample.
3. Find the exterior DMUs with respect to the hull of the best practices.
4. Identify the set of all efficient DMUs.
5. Calculate performance scores as in the traditional DEA model.

This example computes the Big Data radial input-oriented DEA model under variable returns to scale, using random data drawn from a uniform distribution. 500 DMUs with six inputs and four outputs in the interval (10, 20) are generated:
```@example radialbigdata
# Generate random data
using DataEnvelopmentAnalysis
using Distributions
using Random
using StableRNGs

rng = StableRNG(1234567)
X = rand(Uniform(10, 20), 500, 6);
Y = rand(Uniform(10, 20), 500, 4);

# Calculate the Big Data DEA Model
deabig = deabigdata(X, Y)

# Get efficiency scores
efficiency(deabig)
```

### deabigdata Function Documentation

```@docs
deabigdata
```

```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Modified Directional Distance Function

Based on the data  matrix $(X,Y)$, we calculate the modified directional distance function **MDDF**, (Aparicio et al. 2013), of each observation *o* by solving $n$ times the following linear programming problem:

```math
\begin{aligned}
  & \underset{\beta^x, \beta^y,\lambda_j }{\mathop{\min }}\,\quad \quad \quad \;\ \beta^x + \beta^y   \\
  & \text{subject}\ \text{to}  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j x_{ij} }\ \le x_{io} - \beta^x {g}_{io}^{-} \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j y_{rj} }\ \ge y_{ro} + \beta^y {g}_{ro}^{+} \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \lambda_j \ge 0 \qquad j = 1,...,n \\ 
  & \quad \quad \quad \quad \quad \ \beta^x \ge 0 \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \beta^y \ge 0 \qquad r = 1,...,s.  \\
\end{aligned}
```

with the following condition when assuming variable returns to scale:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the modified directional distance function model under variable returns to scale using ones as directions for both inputs and outputs::
```@example mddf
using DataEnvelopmentAnalysis

X = [2; 4; 8; 12; 6; 14; 14; 9.412];

Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example mddf
deamddfvrs= deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS);
nothing # hide
```

```@example mddf
efficiency(deamddfvrs)
```

Estimated $\beta$ on inputs and outputs are returned with the `efficiency` function:
```@example mddf
efficiency(deamddfvrs, :X)
```

```@example mddf
efficiency(deamddfvrs, :Y)
```

### deamddf Function Documentation

```@docs
deamddf
```

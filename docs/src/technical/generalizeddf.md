```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Generalized Distance Function Models

*Chavas and Cox (1999)* introduced a generalized distance function efficiency measure that reescales both inputs and outputs toward the frontier technology.

```math
\begin{aligned}
 & \underset{\delta ,\mathbf{\lambda }}{\mathop{\min }}\,\quad \quad \quad \;\ \delta  \\
 & \text{subject}\ \text{to} \\
 & \quad \quad \quad \quad \quad \ X\mathbf{\lambda } \le \delta^{1 - \alpha} {{\mathbf{x}}_{o}} \\
 & \quad \quad \quad \quad \quad  \;Y\mathbf{\lambda }\ \ge {{\mathbf{y}}_{o}} / \delta^{\alpha} \\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. 
\end{aligned}
```

The measurement of technical efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the generalized distance function DEA model under variable returns to scale using $0.5$ for the value of $\alpha$:
```@example gdf
using DataEnvelopmentAnalysis

X = [5 3; 2 4; 4 2; 4 8; 7 9];

Y = [7 4; 10 8; 8 10; 5 4; 3 6];

deagdf(X, Y, alpha = 0.5, rts = :VRS, slack = false)
```

### deagdf Function Documentation

```@docs
deagdf
```

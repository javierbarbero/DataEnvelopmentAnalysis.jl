```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
    # Solve nonlinear problem to display Ipopt initial message
    X = [1; 2; 3];
    Y = [1; 1; 1];
    deagdf(X, Y, 0.5, rts = :VRS)
end
```

# Generalized Distance Function Models

## Generalized Distance Function Model

*Chavas and Cox (1999)* introduced a generalized distance function efficiency measure that reescales both inputs and outputs toward the frontier technology.

```math
\begin{align}
\label{eq:rim}
  & \underset{\delta ,\mathbf{\lambda }}{\mathop{\min }}\,\quad \quad \quad \;\ \delta  \\
 & \text{subject}\ \text{to} \nonumber \\
 & \quad \quad \quad \quad \quad \ X\mathbf{\lambda } \le \delta^{1 - \alpha} {{\mathbf{x}}_{o}} \nonumber \\
 & \quad \quad \quad \quad \quad  \;Y\mathbf{\lambda }\ \ge {{\mathbf{y}}_{o}} / \delta^{\alpha} \nonumber\\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. \nonumber
\end{align}
```

The measurement of technical efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the generalized distance function DEA model under variable returns to scale using $0.5$ for the value of $\alpha$:
```jldoctest 1
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9];

julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6];

julia> deagdf(X, Y, 0.5, rts = :VRS)
Generalized DF DEA Model 
DMUs = 5; Inputs = 2; Outputs = 2
alpha = 0.5; Returns to Scale = VRS
─────────────────────────────────────────────────────────────────
   efficiency      slackX1      slackX2      slackY1      slackY2
─────────────────────────────────────────────────────────────────
1     0.68185   0.605935     0.0          0.0          4.67865   
2     1.0       0.0          0.0         -1.89292e-7  -1.51434e-7
3     1.0      -3.75555e-8  -1.87777e-8  -7.5111e-8   -9.38887e-8
4     0.25      0.0          0.0         -1.04638e-7  -8.37101e-8
5     0.36      0.2          3.4          3.0         -1.91881e-8
─────────────────────────────────────────────────────────────────
```

### deagdf Function Documentation

```@docs
deagdf
```

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

julia> deagdf(X, Y, 0.5, rts = :VRS, slack = false)
Generalized DF DEA Model
DMUs = 5; Inputs = 2; Outputs = 2
alpha = 0.5; Returns to Scale = VRS
─────────────
   efficiency
─────────────
1     0.68185
2     1.0    
3     1.0    
4     0.25   
5     0.36   
─────────────
```

### deagdf Function Documentation

```@docs
deagdf
```

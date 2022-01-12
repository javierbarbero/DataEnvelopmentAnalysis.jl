```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Radial Models

## Radial Input Oriented Model

Based on the data  matrix $(X,Y)$, we calculate the input oriented efficiency of each observation *o* by solving $n$ times the following linear programming problem -- known as the Charnes, Cooper, and Rhodes (1978), **CCR**, model:
```math
\begin{aligned}
  & \underset{\theta ,\mathbf{\lambda }}{\mathop{\min }}\,\quad \quad \quad \;\ \theta  \\
  & \text{subject}\ \text{to}  \\
  & \quad \quad \quad \quad \quad \ X\mathbf{\lambda } \le \theta {{\mathbf{x}}_{o}} \\
  & \quad \quad \quad \quad \quad  \;Y\mathbf{\lambda }\ \ge {{\mathbf{y}}_{o}}  \\
  & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. 
\end{aligned}
```

The measurement of technical efficiency assuming variable returns to scale, **VRS**, as introduced by *Banker, Charnes and Cooper (1984)* -- known as the Banker, Charnes and Cooper, **BCC**, model -- adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the radial input oriented DEA model under constant returns to scale:
```@example radial
using DataEnvelopmentAnalysis

X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

dea(X, Y, orient = :Input, rts = :CRS)
```

To compute the variable returns to scale model, we simply set the `rts` parameter to `:VRS`:
```@example radial
dea(X, Y, orient = :Input, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example radial
deaiovrs = dea(X, Y, orient = :Input, rts = :VRS);
nothing # hide
```

```@example radial
efficiency(deaiovrs)
```

The optimal peers, ``λ``, are returned with the `peers` function and are returned as a `DEAPeers` object:
```@example radial
peers(deaiovrs)
```

Input and output slacks are returned with the `slacks` function:
```@example radial
slacks(deaiovrs, :X)
```

```@example radial
slacks(deaiovrs, :Y)
```

Input and output optimal targets are returned with the `targets` function:
```@example radial
targets(deaiovrs, :X)
```

```@example radial
targets(deaiovrs, :Y)
```

## Radial Output Oriented Model

It is possible to calculate the output oriented technical efficiency of each observation by solving the following linear program:
```math
\begin{aligned}
 & \underset{\phi ,\mathbf{\lambda }}{\mathop{\max }}\,\quad \quad \quad \quad \phi  \\
 & \text{subject}\ \text{to} \\
 & \quad \quad \quad \quad \quad \ X\lambda\le {{\mathbf{x}}_{o}} \\
 & \quad \quad \quad \quad \quad \ Y\mathbf{\lambda }\ \ge \phi {{\mathbf{y}}_{o}} \\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}.\  & \quad 
\end{aligned}
```

with the following condition when assuming variable returns to scale:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```
In this example we compute the radial output oriented DEA model under variable returns to scale:
```@example radial
dea(X, Y, orient = :Output, rts = :VRS)
```

## Radial Model in Multiplier Form

The dual to the input oriented and output oriented radial DEA models in envelopment form presented above is the multiplier form.

This example computes the radial input-oriented DEA model in multiplier form under constant returns to scale:
```@example radial
deaiovrsm = deam(X, Y, rts = :VRS)
```

Input and output virtual multipliers (shadow prices) are returned with the `multipliers` function:
```@example radial
multipliers(deaiovrsm, :X)
```

```@example radial
multipliers(deaiovrsm, :Y)
```

The value measuring the returns to scale is returned with the `rts` function:
```@example radial
rts(deaiovrsm)
```

### dea Function Documentation

```@docs
dea
deam
```



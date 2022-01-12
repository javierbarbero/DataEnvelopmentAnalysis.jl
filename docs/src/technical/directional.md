```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Directional Distance Function Models

*Chambers, Chung and Fare (1996)* introduced a measure of efficiency that projects observation $\left( {{\mathbf{x}_o,\mathbf{y}_{o}}} \right)$
in a pre-assigned  direction  $\mathbf{g}= {\left({-{\mathbf{g_{x}^-},\mathbf{g^{+}_y}}} \right)\neq\mathbf{0}_{m+s}}$, $\mathbf{g^{-}_{x}}\mathbb{\in R}^m$ and  $\mathbf{g^{+}_{y}}\mathbb{\in R}^s$, in a proportion $\beta$. The associated linear program is:

```math
\begin{aligned}
 & \underset{\beta ,\mathbf{\lambda }}{\mathop{\max }}\,\quad \quad \quad \quad \beta  \\
 & \text{subject}\ \text{to} \\
 & \quad \quad \quad \quad \quad \ X\lambda\le {{\mathbf{x}}_{o}} -\beta{{\mathbf{g^-_x}}} \\
 & \quad \quad \quad \quad \quad \  Y\mathbf{\lambda }\ge\ {{\mathbf{y}}_{o}}+\beta {{\mathbf{g^+_y}}}  \\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}.\\  & \quad 
\end{aligned}
```

The measurement of technical efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the directional distance function DEA model under constant returns to scale using ones as directions for both inputs and outputs:
```@example ddf
using DataEnvelopmentAnalysis

X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

deaddf(X, Y, Gx = :Ones, Gy = :Ones)
```

To compute the variable returns to scale model, we simply set the `rts` parameter to `:VRS`:
```@example ddf
deaddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example ddf
deaddfvrs = deaddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS);
nothing # hide
```

```@example ddf
efficiency(deaddfvrs)
```

The optimal peers, ``λ``, are returned with the `peers` function and are returned as a `DEAPeers` object:
```@example ddf
peers(deaddfvrs)
```
### deaddf Function Documentation

```@docs
deaddf
```

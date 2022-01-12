```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Additive Models

## Weighted Additive Model

The additive model measures technical efficiency based solely on input excesses and output shortfalls, and characterizes efficiency in terms of the input and output slacks: ``\mathbf{s}^-\mathbb{\in R}^m`` and ``\mathbf{s}^+$$\mathbb{\in R}^s``, respectively.
. The package implements the weighted additive formulation of *Cooper and Pastor (1995)* and *Pastor, Lovell and Aparicio (2011)*, whose associated linear program is:
```math
\begin{aligned}
  & \underset{\mathbf{\lambda },\,{{\mathbf{s}}^{-}},\,{{\mathbf{s}}^{+}}}{\mathop{\max }}\,\quad \quad \quad \quad \omega =\mathbf{\rho_{x}^{-}}{{\mathbf{s}}^{\mathbf{-}}}+\mathbf{\rho_{y}^{+}}{{\mathbf{s}}^{+}} \\
 & \text{subject}\ \text{to} \\
 & \quad \quad \quad \quad \quad \quad X\mathbf{\lambda }+{{\mathbf{s}}^{\mathbf{-}}}= \ {{\mathbf{x}}_{o}} \\
 & \quad \quad \quad \quad \quad \quad Y\mathbf{\lambda }-{{\mathbf{s}}^{+}}=\ {{\mathbf{y}}_{o}} \\
 & \quad \quad \quad \quad \quad \quad \mathbf{e\lambda=1} \\
 & \quad \quad \quad \quad \quad \quad \mathbf{\lambda }\ge \mathbf{0},\ {{\mathbf{s}}^{\mathbf{-}}}\ge 0,{{\mathbf{s}}^{+}}\ge 0, 
\end{aligned}
```

where ``(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}}) \mathbb{\in R}^m_{+}\times \mathbb{R}_+^{s}`` are the inputs and outputs weight vectors whose elements can vary across DMUs.

In this example we compute the additive DEA model with all weights equal to one:
```@example additive
using DataEnvelopmentAnalysis

X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

deaadd(X, Y)
```

The same model is computed with:
```@example additive
deaadd(X, Y, :Ones)
```

The additive DEA model can be computed under constant returns to scale setting the
`rts` parameter to `:CRS`:
```@example additive
deaadd(X, Y, :Ones, rts = :CRS)
```

The package can compute a wide class of different DEA models known as general  efficiency measures (GEMs):
- The measure of inefficiency proportions (MIP).
- The normalized weighted additive DEA model.
- The range adjusted measure (RAM).
- The bounded adjusted  measure (BAM).

## Measure of Inefficiency Proportions (MIP)

The measure of inefficiency proportions (MIP), *Charnes et al. (1987)* and *Cooper et al. (1999)*, use the weights:
```math
(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}})=(1/{\mathbf{x}}_{o},1/{{\mathbf{y}}_{o}})
```
```@example additive
deaadd(X, Y, :MIP)
```

## Normalized Weighted Additive Model

The normalized weighted additive DEA model, *Lovell and Pastor (1995)*, use the weights:
```math
(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}})=(1/{\mathbf{σ^-}},1/{{\mathbf{σ^+}}})
```
where ``\mathbf{σ^-}``and ``\mathbf{σ^+}`` are the standard deviations of inputs and outputs respectively.
```@example additive
deaadd(X, Y, :Normalized)
```

## Range Adjusted Measure (RAM)

The range adjusted measure (RAM), *Cooper et al. (1999)*, use the weights::
```math
(\mathbf{\rho^{-}, \mathbf{\rho^{+}}})=(1/(m+s)R^-,(1/(m+s)R^+)
```
where ``R^-``and ``R^+``are the inputs and outputs variables' ranges.
```@example additive
deaadd(X, Y, :RAM)
```

## Bounded Adjusted  Measure (BAM)
The bounded adjusted  measure (BAM), *Cooper et al. (2011)*, use the weights:::
```math
(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}})=(1/(m+s)({\mathbf{x}}_{o}-{\mathbf{\underline{x}}}),(1/(m+s)({\mathbf{\overline{y}}} - {\mathbf{y}}_{o})
```
where ``\mathbf{\underline{x}}`` and ``\mathbf{\overline{y}}`` are the minimum and maximum observed values of inputs and outputs respectively.
```@example additive
deaadd(X, Y, :BAM)
```

### deaadd Function Documentation

```@docs
deaadd
```

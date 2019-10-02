```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Additive Models

## Weighted Additive Model

The additive model measures technical efficiency based solely on input excesses and output shortfalls, and characterizes efficiency in terms of the input and output slacks: ``\mathbf{s}^-\mathbb{\in R}^m`` and ``\mathbf{s}^+$$\mathbb{\in R}^s``, respectively.
. The package implements the weighted additive formulation of *Cooper and Pastor (1995)* and *Pastor, Lovell and Aparicio (2011)*, whose associated linear program is:
```math
\begin{align}
\label{eq:add}
  & \underset{\mathbf{\lambda },\,{{\mathbf{s}}^{-}},\,{{\mathbf{s}}^{+}}}{\mathop{\max }}\,\quad \quad \quad \quad \omega =\mathbf{\rho_{x}^{-}}{{\mathbf{s}}^{\mathbf{-}}}+\mathbf{\rho_{y}^{+}}{{\mathbf{s}}^{+}} \\
 & \text{subject}\ \text{to} \nonumber\\
 & \quad \quad \quad \quad \quad \quad X\mathbf{\lambda }+{{\mathbf{s}}^{\mathbf{-}}}= \ {{\mathbf{x}}_{o}} \nonumber\\
 & \quad \quad \quad \quad \quad \quad Y\mathbf{\lambda }-{{\mathbf{s}}^{+}}=\ {{\mathbf{y}}_{o}} \nonumber\\
 & \quad \quad \quad \quad \quad \quad \mathbf{e\lambda=1} \nonumber\\
 & \quad \quad \quad \quad \quad \quad \mathbf{\lambda }\ge \mathbf{0},\ {{\mathbf{s}}^{\mathbf{-}}}\ge 0,{{\mathbf{s}}^{+}}\ge 0, \nonumber
\end{align}
```

where ``(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}}) \mathbb{\in R}^m_{+}\times \mathbb{R}_+^{s}`` are the inputs and outputs weight vectors whose elements can vary across DMUs.

In this example we compute the additive DEA model with all weights equal to one:
```jldoctest 1
julia> using DataEnvelopmentAnalysis

julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> deaadd(X, Y)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = Ones; Returns to Scale = VRS
────────────────────────────────────────────────────
      efficiency       slackX1  slackX2      slackY1
────────────────────────────────────────────────────
1    0.0           0.0              0.0  0.0
2    7.33333       4.33333          0.0  3.0
3    0.0           0.0              0.0  0.0
4   -8.03397e-16  -8.03397e-16      0.0  0.0
5   18.0          13.0              1.0  4.0
6    6.48305e-16   2.70127e-16      0.0  3.78178e-16
7    0.0           0.0              0.0  0.0
8    0.0           0.0              0.0  0.0
9    0.0           0.0              0.0  0.0
10  35.0          25.0             10.0  0.0
11   4.0           0.0              4.0  4.78849e-16
────────────────────────────────────────────────────
```

The same model is computed with:
```julia
deaadd(X, Y, :Ones)
```

The additive DEA model can be computed under constant returns to scale setting the
`rts` parameter to `:CRS`:
```julia
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
```jldoctest 1
julia> deaadd(X, Y, :MIP)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = MIP; Returns to Scale = VRS
─────────────────────────────────────────────────────
      efficiency       slackX1  slackX2       slackY1
─────────────────────────────────────────────────────
1    0.0           0.0              0.0   0.0
2    0.507519      0.0              0.0   7.10526
3    0.0           0.0              0.0   0.0
4   -4.72586e-17  -8.03397e-16      0.0   0.0
5    2.20395       0.0              0.0  17.6316
6    1.31279e-16   8.10382e-16      0.0   8.64407e-16
7    0.0           0.0              0.0   0.0
8    0.0           0.0              0.0   0.0
9    0.0           0.0              0.0   0.0
10   1.04322      17.0             15.0   1.0
11   0.235294      0.0              4.0   0.0
─────────────────────────────────────────────────────
```

## Normalized Weighted Additive Model

The normalized weighted additive DEA model, *Lovell and Pastor (1995)*, use the weights:
```math
(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}})=(1/{\mathbf{σ^-}},1/{{\mathbf{σ^+}}})
```
where ``\mathbf{σ^-}``and ``\mathbf{σ^+}`` are the standard deviations of inputs and outputs respectively.
```jldoctest 1
julia> deaadd(X, Y, :Normalized)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = Normalized; Returns to Scale = VRS
──────────────────────────────────────────────────────────
      efficiency       slackX1       slackX2       slackY1
──────────────────────────────────────────────────────────
1    0.0           0.0           0.0           0.0
2    0.804925      0.0           0.65          6.25
3    0.0           0.0           0.0           0.0
4   -9.79609e-17   0.0          -6.09909e-16   0.0
5    2.01497       0.0           2.95         13.75
6    4.81529e-16   2.49057e-15   0.0           2.37658e-15
7    0.0           0.0           0.0           0.0
8    0.0           0.0           0.0           0.0
9    0.0           0.0           0.0           0.0
10   3.98989      17.0          15.0           1.0
11   0.642462      0.0           4.0           0.0
──────────────────────────────────────────────────────────
```

## Range Adjusted Measure (RAM)

The range adjusted measure (RAM), *Cooper et al. (1999)*, use the weights::
```math
(\mathbf{\rho^{-}, \mathbf{\rho^{+}}})=(1/(m+s)R^-,(1/(m+s)R^+)
```
where ``R^-``and ``R^+``are the inputs and outputs variables' ranges.
```jldoctest 1
julia> deaadd(X, Y, :RAM)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = RAM; Returns to Scale = VRS
──────────────────────────────────────────────────────────
      efficiency       slackX1       slackX2       slackY1
──────────────────────────────────────────────────────────
1    0.0           0.0           0.0           0.0
2    0.102975      0.0           0.0           7.10526
3    0.0           0.0           0.0           0.0
4   -1.01651e-17   0.0          -6.09909e-16   0.0
5    0.25553       0.0           0.0          17.6316
6    5.68808e-17   2.49057e-15   0.0           2.37658e-15
7    0.0           0.0           0.0           0.0
8    0.0           0.0           0.0           0.0
9    0.0           0.0           0.0           0.0
10   0.417646     17.0          15.0           1.0
11   0.0666667     0.0           4.0           0.0
──────────────────────────────────────────────────────────
```

## Bounded Adjusted  Measure (BAM)
The bounded adjusted  measure (BAM), *Cooper et al. (2011)*, use the weights:::
```math
(\mathbf{\rho_{x}^{-}, \mathbf{\rho_{y}^{+}}})=(1/(m+s)({\mathbf{x}}_{o}-{\mathbf{\underline{x}}}),(1/(m+s)({\mathbf{\overline{y}}} - {\mathbf{y}}_{o})
```
where ``\mathbf{\underline{x}}`` and ``\mathbf{\overline{y}}`` are the minimum and maximum observed values of inputs and outputs respectively.
```jldoctest 1
julia> deaadd(X, Y, :BAM)
Weighted Additive DEA Model
DMUs = 11; Inputs = 2; Outputs = 1
Weights = BAM; Returns to Scale = VRS
─────────────────────────────────────────────────
      efficiency   slackX1  slackX2       slackY1
─────────────────────────────────────────────────
1    0.0           0.0          0.0   0.0
2    0.199894      6.59649      0.0   0.0
3    0.0           0.0          0.0   0.0
4   -3.78838e-17   0.0          0.0  -5.68256e-16
5    0.432971     13.0          1.0   4.0
6    1.11554e-17   0.0          0.0   7.36254e-16
7    0.0           0.0          0.0   0.0
8    0.0           0.0          0.0   0.0
9    0.0           0.0          0.0   0.0
10   0.571361      5.0         11.0   5.0
11   0.121212      0.0          4.0   0.0
─────────────────────────────────────────────────
```

### deaadd Function Documentation

```@docs
deaadd
```

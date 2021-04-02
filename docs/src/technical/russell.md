```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
    # Solve nonlinear problem to display Ipopt initial message
    dearussell([1; 2; 3], [1; 1; 1], orient = :Graph, rts = :VRS)
end
```

# Russell Models

## Russell Input Model

Based on the data  matrix $(X,Y)$, we calculate the Russell measure of input efficiency (Färe & Lovell, 1978; and Färe et al., 1985) of each observation *o* by solving $n$ times the following linear programming problem:
```math
\begin{aligned}
  & \underset{\theta_i ,\lambda_j }{\mathop{\min }}\,\quad \quad \quad \;\ \frac{1}{m} \sum_{i=1}^{m}{\theta_i}  \\
  & \text{subject}\ \text{to}  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j x_{ij} }\ \le \theta_i {x}_{io} \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j y_{rj} }\ \ge {y}_{ro} \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \theta_i \le 1 \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \lambda_j \ge 0 \qquad j = 1,...,n. 
\end{aligned}
```

The measurement of technical efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the Russell input DEA model under constant returns to scale:
```jldoctest 1
julia> using DataEnvelopmentAnalysis

julia> X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6 8];

julia> Y = [1; 1; 1; 1; 1; 1; 1; 1];

julia> dearussell(X, Y, orient = :Input, rts = :CRS)
Russell DEA Model 
DMUs = 8; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = CRS
──────────────────────────────────────────
   efficiency     effX1     effX2  slackY1
──────────────────────────────────────────
1    1.0       1.0       1.0           0.0
2    1.0       1.0       1.0           0.0
3    1.0       1.0       1.0           0.0
4    0.583333  0.5       0.666667      0.0
5    0.4       0.4       0.4           0.0
6    0.833333  0.666667  1.0           0.0
7    0.65      0.5       0.8           0.0
8    0.5625    0.625     0.5           0.0
──────────────────────────────────────────
```

To compute the variable returns to scale model, we simply set the `rts` parameter to `:VRS`:

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> dearussellio = dearussell(X, Y, orient = :Input, rts = :CRS);

julia> efficiency(dearussellio)
8-element Vector{Float64}:
 1.0
 1.0
 1.0
 0.5833333333333334
 0.4
 0.8333333333333334
 0.65
 0.5625

julia> efficiency(dearussellio, :X)
8×2 Matrix{Float64}:
 1.0       1.0
 1.0       1.0
 1.0       1.0
 0.5       0.666667
 0.4       0.4
 0.666667  1.0
 0.5       0.8
 0.625     0.5
```

## Russell Output Model

It is possible to calculate the Russell measure of output efficiency of each observation by solving the following linear program:

```math
\begin{aligned}
  & \underset{\phi_r ,\lambda_j }{\mathop{\max }}\,\quad \quad \quad \;\ \frac{1}{s} \sum_{r=1}^{s}{\phi_r}  \\
  & \text{subject}\ \text{to}  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j x_{ij} }\ \le {x}_{io} \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j y_{rj} }\ \ge \phi_r {y}_{ro} \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \phi_r \ge 1 \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \lambda_j \ge 0 \qquad j = 1,...,n. 
\end{aligned}
```

with the following condition when assuming variable returns to scale:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```
In this example we compute the Russell output DEA model under constant returns to scale:
```jldoctest 1
julia> X = [1; 1; 1; 1; 1; 1; 1; 1];

julia> Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5] ;

julia> dearussell(X, Y, orient = :Output, rts = :CRS)
Russell DEA Model 
DMUs = 8; Inputs = 1; Outputs = 2
Orientation = Output; Returns to Scale = CRS
────────────────────────────────────────
   efficiency    effY1    effY2  slackX1
────────────────────────────────────────
1     1.0      1.0      1.0          0.0
2     1.0      1.0      1.0          0.0
3     1.0      1.0      1.0          0.0
4     1.86667  2.33333  1.4          0.0
5     2.33333  2.33333  2.33333      0.0
6     1.5      1.0      2.0          0.0
7     1.45833  1.16667  1.75         0.0
8     3.05556  5.11111  1.0          0.0
────────────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> dearusselloo = dearussell(X, Y, orient = :Output, rts = :CRS);

julia> efficiency(dearusselloo)
8-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.8666666666666665
 2.333333333333333
 1.5
 1.4583333333333335
 3.0555555555555554

julia> efficiency(dearusselloo, :Y)
8×2 Matrix{Float64}:
 1.0      1.0
 1.0      1.0
 1.0      1.0
 2.33333  1.4
 2.33333  2.33333
 1.0      2.0
 1.16667  1.75
 5.11111  1.0
```

## Russell Graph Model

It is possible to calculate the Russell graph measure of technical efficiency of each observation by solving the following linear program:
```math
\begin{aligned}
  & \underset{\theta_i, \phi_r ,\lambda_j }{\mathop{\min }}\,\quad \quad \quad \;\ \frac{1}{m + s} (\sum_{i=1}^{m}{\theta_i} +  \sum_{r=1}^{s}{\frac{1}{\phi_r}})  \\
  & \text{subject}\ \text{to}  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j x_{ij} }\ \le \theta_i {x}_{io}  \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\lambda_j y_{rj} }\ \ge \phi_r {y}_{ro} \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \theta_i \le 1 \qquad i = 1,...,m  \\
  & \quad \quad \quad \quad \quad \ \phi_r \ge 1 \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \lambda_j \ge 0 \qquad j = 1,...,n. 
\end{aligned}
```

with the following condition when assuming variable returns to scale:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```
In this example we compute the Russell graph DEA model under variable returns to scale:
```jldoctest 1
julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> dearussell(X, Y, orient = :Graph, rts = :VRS)
Russell DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Orientation = Graph; Returns to Scale = VRS
────────────────────────────────
   efficiency     effX1    effY1
────────────────────────────────
1    1.0       1.0       1.0
2    1.0       1.0       1.0
3    1.0       1.0       1.0
4    1.0       1.0       1.0
5    0.633333  0.666667  1.66667
6    0.723214  0.571429  1.14286
7    0.928571  0.857143  1.0
8    0.447795  0.424989  2.12495
────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> dearussellgr = dearussell(X, Y, orient = :Graph, rts = :VRS);

julia> efficiency(dearussellgr)
8-element Vector{Float64}:
 0.9999996498240314
 0.99999997633837
 0.9999999894257894
 0.9999999706938013
 0.6333333226351023
 0.7232142937165499
 0.9285714048846612
 0.4477946881737884

julia> efficiency(dearussellgr, :X)
8×1 Matrix{Float64}:
 1.0
 0.9999999519440551
 0.9999999739201625
 0.9999999335267051
 0.6666666303008444
 0.57142851550475
 0.8571428024451365
 0.4249893701037514

julia> efficiency(dearussellgr, :Y)
8×1 Matrix{Float64}:
 1.0000007087261955
 1.0
 1.0
 1.0
 1.6666666250851112
 1.1428570489099183
 1.0
 2.124946848134728
```

### dearussell Function Documentation

```@docs
dearussell
```

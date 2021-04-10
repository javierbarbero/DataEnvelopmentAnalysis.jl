```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
    # Solve nonlinear problem to display Ipopt initial message
    deamddf([1; 2; 3], [1; 1; 1],Gx = :Ones, Gy = :Ones, rts = :VRS, slack = false)
end
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
```jldoctest 1
julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS)
Modified DDF DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Returns to Scale = VRS
Gx = Ones; Gy = Ones
───────────────────────────────────────────────────────
   efficiency          βx          βy  slackX1  slackY1
───────────────────────────────────────────────────────
1  7.42338e-7  0.0         7.49867e-7      0.0      0.0
2  3.97731e-7  3.78622e-7  1.91088e-8      0.0      0.0
3  4.83998e-7  4.86874e-7  0.0             0.0      0.0
4  1.32755e-6  1.33672e-6  0.0             0.0      0.0
5  4.0         2.0         2.0             0.0      0.0
6  7.33333     7.33333     0.0             0.0      0.0
7  2.0         2.0         0.0             0.0      0.0
8  8.059       5.412       2.647           0.0      0.0
───────────────────────────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> deamddfvrs= deamddf(X, Y, Gx = :Ones, Gy = :Ones, rts = :VRS);

julia> efficiency(deamddfvrs)
8-element Vector{Float64}:
 7.423380868494941e-7
 3.9773060690474414e-7
 4.839979925116349e-7
 1.327547239652975e-6
 4.000000398014332
 7.333333863503837
 2.0000013477402137
 8.059000425674004
```

Estimated $\beta$ on inputs and outputs are returned with the `efficiency` function:
```jldoctest 1
julia> efficiency(deamddfvrs, :X)
8-element Vector{Float64}:
 0.0
 3.786218512295063e-7
 4.868738237228173e-7
 1.3367161141172963e-6
 2.0000004088097754
 7.333333865967325
 2.0000013569193515
 5.412000442949787

julia> efficiency(deamddfvrs, :Y)
8-element Vector{Float64}:
 7.498667608019586e-7
 1.910875567523785e-8
 0.0
 0.0
 1.9999999892045561
 0.0
 0.0
 2.6469999827242163
```

### deamddf Function Documentation

```@docs
deamddf
```

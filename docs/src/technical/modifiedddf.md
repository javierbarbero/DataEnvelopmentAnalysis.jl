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
1  7.42311e-7  0.0         7.49834e-7      0.0      0.0
2  4.03943e-7  3.95267e-7  8.67603e-9      0.0      0.0
3  4.84094e-7  4.87082e-7  0.0             0.0      0.0
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
8-element Array{Float64,1}:
 7.423112595976122e-7
 4.0394299500299253e-7
 4.840938492567068e-7
 1.3275492427779482e-6
 4.000000398014331
 7.333333863503834
 2.0000013699988735
 8.059000425674004
```

Estimated $\beta$ on inputs and outputs are returned with the `efficiency` function:
```jldoctest 1
julia> efficiency(deamddfvrs, :X)
8-element Array{Float64,1}:
 0.0
 3.95266969245659e-7
 4.870818366962988e-7
 1.3367167143310045e-6
 2.000000408809775
 7.333333865967322
 2.0000013799988317
 5.412000442949787

julia> efficiency(deamddfvrs, :Y)
8-element Array{Float64,1}:
 7.498338684288291e-7
 8.67602575733351e-9
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

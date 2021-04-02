```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Enhanced Russell Graph Slack Based Measure

Based on the data  matrix $(X,Y)$, we calculate the Enhanced Russell Graph Measure, **ERG**, (Pastor et al., 1999) -- also known as the Slack Based Measure, **SBM**, Tone (2001) -- of each observation *o* by solving $n$ times the following linear programming problem:

```math
\begin{aligned}
  & \underset{\beta, t_i^-, t_r^+ ,\mu_j }{\mathop{\min }}\,\quad \quad \quad \;\ \beta -  \frac{1}{m} \sum_{i=1}^{m}{\frac{t_i^-}{x_{io}}}  \\
  & \text{subject}\ \text{to}  \\
  & \quad \quad \quad \quad \quad \ \beta +  \frac{1}{s} \sum_{r=1}^{s}{\frac{t_r^+}{y_{ro}}} = 1 \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\mu_j x_{ij} }\ = \beta {x}_{io} - t_i^- \qquad i = 1,...,m \\
  & \quad \quad \quad \quad \quad \ \sum_{j=1}^{n}{\mu_j y_{rj} }\ = \beta {y}_{ro} + t_r^+ \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \beta \ge 0  \\
  & \quad \quad \quad \quad \quad \ t_i^- \ge 0 \qquad i = 1,...,m \\
  & \quad \quad \quad \quad \quad \ t_r^+ \ge 0 \qquad r = 1,...,s \\
  & \quad \quad \quad \quad \quad \ \mu_j \ge 0 \qquad j = 1,...,n. 
\end{aligned}
```

with the following condition when assuming variable returns to scale:
```math
\sum\nolimits_{j=1}^{n}\mu_j=\beta
```

After solving the model, input and output slacks are recovered through the following expressions:
```math
\begin{aligned}
 & s_i^- = \frac{t_i^-}{\beta} \qquad i = 1,...,m \\
 & s_r^+ = \frac{t_r^+}{\beta} \qquad r = 1,...,s.
\end{aligned}
```

In this example we compute the Enhanced Russell Graph DEA model under variable returns to scale:
```jldoctest 1
julia> X = [2; 4; 8; 12; 6; 14; 14; 9.412];

julia> Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

julia> deaerg(X, Y, rts = :VRS)
Enhanced Russell Graph Slack Based Measure DEA Model 
DMUs = 8; Inputs = 1; Outputs = 1
Orientation = Graph; Returns to Scale = VRS
───────────────────────────────────────
   efficiency    beta  slackX1  slackY1
───────────────────────────────────────
1    1.0       1.0     0.0        0.0
2    1.0       1.0     0.0        0.0
3    1.0       1.0     0.0        0.0
4    1.0       1.0     0.0        0.0
5    0.4       0.6     2.0        2.0
6    0.47619   1.0     7.33333    0.0
7    0.857143  1.0     2.0        0.0
8    0.2       0.4706  5.412      2.647
───────────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> deaergvrs = deaerg(X, Y, rts = :VRS);

julia> efficiency(deaergvrs)
8-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0
 0.39999999999999997
 0.47619047619047616
 0.8571428571428574
 0.2
```

Estimated $\beta$'s are returned with the `efficiency` function using `:beta` as the second argument:
```jldoctest 1
julia> efficiency(deaergvrs, :beta)
8-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0
 0.6
 0.9999999999999998
 0.9999999999999998
 0.4706000000000001
```

### deaerg Function Documentation

```@docs
deaerg
```

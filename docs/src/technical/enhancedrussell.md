```@meta
CurrentModule = DataEnvelopmentAnalysis
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
```@example ergsbm
using DataEnvelopmentAnalysis

X = [2; 4; 8; 12; 6; 14; 14; 9.412];

Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

deaergvrs = deaerg(X, Y, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example ergsbm
efficiency(deaergvrs)
```

Estimated $\beta$'s are returned with the `efficiency` function using `:beta` as the second argument:
```@example ergsbm
efficiency(deaergvrs, :beta)
```

### deaerg Function Documentation

```@docs
deaerg
```

```@meta
CurrentModule = DataEnvelopmentAnalysis
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
```@example russell
using DataEnvelopmentAnalysis

X = [2 2; 1 4; 4 1; 4 3; 5 5; 6 1; 2 5; 1.6 8];

Y = [1; 1; 1; 1; 1; 1; 1; 1];

dearussell(X, Y, orient = :Input, rts = :CRS)
```

To compute the variable returns to scale model, we simply set the `rts` parameter to `:VRS`:

Estimated efficiency scores are returned with the `efficiency` function:
```@example russell
dearussellio = dearussell(X, Y, orient = :Input, rts = :CRS);
nothing # hide
```

```@example russell
efficiency(dearussellio)
```

```@example russell
efficiency(dearussellio, :X)
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
```@example russell
X = [1; 1; 1; 1; 1; 1; 1; 1];

Y = [7 7; 4 8; 8 4; 3 5; 3 3; 8 2; 6 4; 1.5 5] ;

dearussell(X, Y, orient = :Output, rts = :CRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example russell
dearusselloo = dearussell(X, Y, orient = :Output, rts = :CRS);
nothing # hide
```
```@example russell
efficiency(dearusselloo)
```

```@example russell
efficiency(dearusselloo, :Y)
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
```@example russell
X = [2; 4; 8; 12; 6; 14; 14; 9.412];

Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

dearussell(X, Y, orient = :Graph, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example russell
dearussellgr = dearussell(X, Y, orient = :Graph, rts = :VRS);
nothing # hide
```

```@example russell
efficiency(dearussellgr)
```

```@example russell
efficiency(dearussellgr, :X)
```
```@example russell
efficiency(dearussellgr, :Y)
```

### dearussell Function Documentation

```@docs
dearussell
```

```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Radial Models

## Radial Input Oriented Model

Based on the data  matrix $(X,Y)$, we calculate the input oriented efficiency of each observation *o* by solving $n$ times the following linear programming problem -- known as the Charnes, Cooper, and Rhodes (1978), **CCR**, model:
```math
\begin{align}
\label{eq:rim}
  & \underset{\theta ,\mathbf{\lambda }}{\mathop{\min }}\,\quad \quad \quad \;\ \theta  \\
 & \text{subject}\ \text{to} \nonumber \\
 & \quad \quad \quad \quad \quad \ X\mathbf{\lambda } \le \theta {{\mathbf{x}}_{o}} \nonumber \\
 & \quad \quad \quad \quad \quad  \;Y\mathbf{\lambda }\ \ge {{\mathbf{y}}_{o}}  \nonumber\\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. \nonumber
\end{align}
```
The measurement of technical efficiency assuming variable returns to scale, **VRS**, as introduced by *Banker, Charnes and Cooper (1984)* -- known as the Banker, Charnes and Cooper, **BCC**, model -- adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

In this example we compute the radial input oriented DEA model under constant returns to scale:
```jldoctest 1
julia> using DataEnvelopmentAnalysis

julia> X = [5 13; 16 12; 16 26; 17 15; 18 14; 23 6; 25 10; 27 22; 37 14; 42 25; 5 17];

julia> Y = [12; 14; 25; 26; 8; 9; 27; 30; 31; 26; 12];

julia> dea(X, Y, orient = :Input, rts = :CRS)
Radial DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = CRS
──────────────────────────────────────────────────
    efficiency       slackX1      slackX2  slackY1
──────────────────────────────────────────────────
1     1.0        0.0          0.0              0.0
2     0.62229   -4.41868e-15  0.0              0.0
3     0.819856   0.0          8.17926e-15      0.0
4     1.0       -8.03397e-16  0.0              0.0
5     0.310371   1.80764e-15  0.0              0.0
6     0.555556   4.44444      0.0              0.0
7     1.0        0.0          0.0              0.0
8     0.757669   1.60679e-15  0.0              0.0
9     0.820106   1.64021      0.0              0.0
10    0.490566   9.68683e-15  0.0              0.0
11    1.0        0.0          4.0              0.0
──────────────────────────────────────────────────
```

To compute the variable returns to scale model, we simply set the `rts` parameter to `:VRS`:
```jldoctest 1
julia> dea(X, Y, orient = :Input, rts = :VRS)
Radial DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Input; Returns to Scale = VRS
───────────────────────────────────────────────────────
    efficiency       slackX1       slackX2      slackY1
───────────────────────────────────────────────────────
1     1.0        0.0           0.0          0.0        
2     0.869986   0.0           0.0          0.0        
3     1.0        0.0           2.56789e-13  0.0        
4     1.0       -8.03397e-16   0.0          0.0        
5     0.71164    0.0           0.0          2.69841    
6     1.0        2.70127e-16   0.0          3.78178e-16
7     1.0        0.0           0.0          0.0        
8     1.0        0.0          -1.27018e-14  0.0        
9     1.0        0.0           0.0          0.0        
10    0.493121   3.90444e-15   0.0          0.0        
11    1.0        0.0           4.0          4.78849e-16
───────────────────────────────────────────────────────
```

Estimated efficiency scores are returned with the `efficiency` function:
```jldoctest 1
julia> deaiovrs = dea(X, Y, orient = :Input, rts = :VRS);

julia> efficiency(deaiovrs)
11-element Array{Float64,1}:
 1.0
 0.8699861687413553
 1.0000000000000002
 1.0
 0.7116402116402116
 1.0
 1.0
 0.9999999999999999
 1.0
 0.4931209268645909
 1.0
```

The optimal peers, ``λ``, are returned with the `peers` function and are returned as a `SparseArrays.SparseMatrixCSC{Float64,Int64}` object:
```julia
peers(deaiovrs)
```


## Radial Output Oriented Model

It is possible to calculate the output oriented technical efficiency of each observation by solving the following linear program:
```math
\begin{align}
\label{eq:rom}
  & \underset{\phi ,\mathbf{\lambda }}{\mathop{\max }}\,\quad \quad \quad \quad \phi  \\
 & \text{subject}\ \text{to} \nonumber\\
 & \quad \quad \quad \quad \quad \ X\lambda\le {{\mathbf{x}}_{o}} \nonumber\\
 & \quad \quad \quad \quad \quad \ Y\mathbf{\lambda }\ \ge \phi {{\mathbf{y}}_{o}}  \nonumber\\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}.\nonumber\\  & \quad \nonumber
\end{align}
```

with the following condition when assuming variable returns to scale:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```
In this example we compute the radial output oriented DEA model under variable returns to scale:
```jldoctest 1
julia> dea(X, Y, orient = :Output, rts = :VRS)
Radial DEA Model 
DMUs = 11; Inputs = 2; Outputs = 1
Orientation = Output; Returns to Scale = VRS
──────────────────────────────────────────────────
    efficiency       slackX1  slackX2      slackY1
──────────────────────────────────────────────────
1      1.0       0.0              0.0  0.0        
2      1.50752   5.78599e-15      0.0  0.0        
3      1.0       0.0              0.0  0.0        
4      1.0      -8.03397e-16      0.0  0.0        
5      3.20395  -3.38377e-15      0.0  0.0        
6      1.0       2.70127e-16      0.0  3.78178e-16
7      1.0       0.0              0.0  0.0        
8      1.0       0.0              0.0  0.0        
9      1.0       0.0              0.0  0.0        
10     1.19231   5.0             11.0  0.0        
11     1.0       0.0              4.0  4.78849e-16
──────────────────────────────────────────────────
```
### dea Function Documentation

```@docs
dea
```

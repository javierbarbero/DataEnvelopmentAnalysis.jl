```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Profitability Models

## Profitability Model

The profitabilty function defines as $\mathrm{P}\left(\mathbf{w},\mathbf{p}\right)=\max \Big\{ \sum\limits_{i=1}^{s}{{p}_{i}}{{y}_{i}}/\sum\limits_{i=1}^{m}{{w}_{i}}{{x}_{i}} \,|  {\mathbf{x}} \geqslant X\mathbf{\lambda},\;{\mathbf{y}} \leqslant Y{\mathbf{\lambda },\; \lambda } \geqslant {\mathbf{0}} \Big\}$. *Zofío and Prieto (2006)* introduced the following program that allows calculating profitability efficiency.
```math
\begin{align}
\label{eq:maxprofit}
  & \underset{\mathbf{x,y,\lambda_{j},\omega} }{\mathop{\min }}\,\quad \quad \quad \;\ \omega  \\ 
 & \text{subject}\ \text{to} \nonumber \\ 
 & \quad \quad \quad \quad \quad \ {\sum_{j=1}^{j} \lambda^{j} \frac{w^{j} x^{j}}{p^{j} y^{j}} = \omega \frac{w^{j} x^{j}_{o}}{p^{j} y^{j}_{o}} } \nonumber \\ 
 & \quad \quad \quad \quad \quad  \; \sum\nolimits_{j=1}^{n}\lambda^{j}=1 \nonumber\\

\nonumber \\ 
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. \nonumber  
\end{align}
```
*Profitabilty efficiency* defines as the ratio between maximum profitabilty and observed profitabilty. Following the duality results introduced by *Zofío and Prieto (2006)* it is possible to decompose it into technical and allocative efficiencies under constant returns to scale. Profitabilty efficiency can be then decomposed into the generalizaed distance fucntion and the residual ratio corresponding to the *allocative profit efficiency*. Allocative efficiency defines then as the ratio of profitability at the technically efficient projection on the frontier to maximum profitability. 

In this example we compute the profitability efficiency measure:
```jldoctest 1
julia> X = [5 3; 2 4; 4 2; 4 8; 7 9.0];

julia> Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];

julia> W = [2 1; 2 1; 2 1; 2 1; 2 1.0];

julia> P = [3 2; 3 2; 3 2; 3 2; 3 2.0];

julia> deaprofitability(X, Y, W, P)
Profitability DEA Model
DMUs = 5; Inputs = 2; Outputs = 2
alpha = 0.5; Returns to Scale = VRS
─────────────────────────────────────────────────────────
   Profitability       CRS      VRS     Scale  Allocative
─────────────────────────────────────────────────────────
1       0.38796   0.636364  0.68185  0.93329     0.609651
2       1.0       1.0       1.0      1.0         1.0
3       0.765217  1.0       1.0      1.0         0.765217
4       0.25      0.25      0.25     1.0         1.0
5       0.15879   0.26087   0.36     0.724638    0.608696
─────────────────────────────────────────────────────────
```

### deaprofitability Function Documentation

```@docs
deaprofitability
```

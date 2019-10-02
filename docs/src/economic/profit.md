```@meta
CurrentModule = DataEnvelopmentAnalysis
DocTestSetup = quote
    using DataEnvelopmentAnalysis
end
```

# Profit Models

## Profit Efficiency Model with Directional Distance Function Technical Efficiency

 The profit function defines as $\Pi\left(\mathbf{w},\mathbf{p}\right)=\max \Big\{ \sum\limits_{i=1}^{s}{{p}_{i}}{{y}_{i}}-\sum\limits_{i=1}^{m}{{w}_{i}}{{x}_{i}} \,|  $
 $ {\mathbf{x}} \geqslant X\mathbf{\lambda},\;{\mathbf{y}} \leqslant Y{\mathbf{\lambda },\;{\mathbf{\mathbf{e\lambda=1}, \lambda }} \geqslant {\mathbf{0}}} \Big\}$. Calculating maximum profit along with the optimal output and input quantities $\mathbf{y^{*}}$and $\mathbf{x^{*}}$ requires solving: 

```math
\begin{align}
\label{eq:maxprofit}
  & \underset{\mathbf{x,y,\lambda} }{\mathop{\max }}\,\quad \quad \quad \;\ \Pi\left(\mathbf{w},\mathbf{p}\right)=\mathbf{py^{*}-wx^{*}}  \\ 
 & \text{subject}\ \text{to} \nonumber \\ 
 & \quad \quad \quad \quad \quad \ {{\mathbf{x}}}\ge X\mathbf{\lambda=x } \nonumber \\ 
 & \quad \quad \quad \quad \quad  \; {{\mathbf{y}}}  \le Y\mathbf{\lambda =y} \nonumber\\
& \quad \quad \quad \quad \quad \; \mathbf{e\lambda=1}
\nonumber \\ 
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. \nonumber  
\end{align}
```

*Profit efficiency* defines as the difference between maximum profit and observed profit. Following the duality results introduced by *Chambers, Chung and Färe (1998)* it is possible to decompose it into technical and allocative efficiencies under variable returns to scale. Profit efficiency can be then decomposed into the directional distance fucntion and the residual difference corresponding to the *allocative profit efficiency*. Allocative efficiency defines then as the difference between maximum profit and profit at the technically efficient projection on the frontier. The approach relies on the directional vector to normalize these components, thereby ensuring that their values can be compared across DMUs. 

In this example we compute the profit efficiency measure under variable returns to scale:
```jldoctest 1
julia> X = [1 1; 1 1; 0.75 1.5; 0.5 2; 0.5 2; 2 2; 2.75 3.5; 1.375 1.75];

julia> Y = [1 11; 5 3; 5 5; 2 9; 4 5; 4 2; 3 3; 4.5 3.5];

julia> P = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

julia> W = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

julia> GxGydollar = 1 ./ (sum(P, dims = 2) + sum(W, dims = 2));

julia> Gx = repeat(GxGydollar, 1, 2);

julia> Gy = repeat(GxGydollar, 1, 2);

julia> deaprofit(X, Y, W, P, Gx, Gy)
Profit DEA Model 
DMUs = 8; Inputs = 2; Outputs = 2
Returns to Scale = VRS
─────────────────────────────────────
   Profit     Technical    Allocative
─────────────────────────────────────
1     2.0   0.0           2.0        
2     2.0  -5.41234e-16   2.0        
3     0.0   0.0           0.0        
4     2.0   0.0           2.0        
5     2.0   0.0           2.0        
6     8.0   6.0           2.0        
7    12.0  12.0          -1.77636e-15
8     4.0   3.0           1.0        
─────────────────────────────────────
```

### deaprofit Function Documentation

```@docs
deaprofit
```


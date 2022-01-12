```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Profitability Models

## Profitability Model

The profitabilty function defines as $\mathrm{P}\left(\mathbf{w},\mathbf{p}\right)=\max \Big\{ \sum\limits_{i=1}^{s}{{p}_{i}}{{y}_{i}}/\sum\limits_{i=1}^{m}{{w}_{i}}{{x}_{i}} \,|  {\mathbf{x}} \geqslant X\mathbf{\lambda},\;{\mathbf{y}} \leqslant Y{\mathbf{\lambda },\; \lambda } \geqslant {\mathbf{0}} \Big\}$. *Zofío and Prieto (2006)* introduced the following program that allows calculating profitability efficiency.
```math
\begin{aligned}
 & \underset{\mathbf{x,y,\lambda_{j},\omega} }{\mathop{\min }}\,\quad \quad \quad \;\ \omega  \\
 & \text{subject}\ \text{to} \\
 & \quad \quad \quad \quad \quad \ {\sum_{j=1}^{j} \lambda^{j} \frac{w^{j} x^{j}}{p^{j} y^{j}} = \omega \frac{w^{j} x^{j}_{o}}{p^{j} y^{j}_{o}} } \\
 & \quad \quad \quad \quad \quad  \; \sum\nolimits_{j=1}^{n}\lambda^{j}=1 \\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. 
\end{aligned}
```
*Profitabilty efficiency* defines as the ratio between maximum profitabilty and observed profitabilty. Following the duality results introduced by *Zofío and Prieto (2006)* it is possible to decompose it into technical and allocative efficiencies under constant returns to scale. Profitabilty efficiency can be then decomposed into the generalizaed distance fucntion and the residual ratio corresponding to the *allocative profit efficiency*. Allocative efficiency defines then as the ratio of profitability at the technically efficient projection on the frontier to maximum profitability.

In this example we compute the profitability efficiency measure:
```@example profitability
using DataEnvelopmentAnalysis

X = [5 3; 2 4; 4 2; 4 8; 7 9.0];

Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];

W = [2 1; 2 1; 2 1; 2 1; 2 1.0];

P = [3 2; 3 2; 3 2; 3 2; 3 2.0];

deaprofitability(X, Y, W, P)
```

### deaprofitability Function Documentation

```@docs
deaprofitability
```

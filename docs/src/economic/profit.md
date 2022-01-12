```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Profit Models

## Profit Efficiency Model with Directional Distance Function Technical Efficiency

 The profit function defines as $\Pi\left(\mathbf{w},\mathbf{p}\right)=\max \Big\{ \sum\limits_{i=1}^{s}{{p}_{i}}{{y}_{i}}-\sum\limits_{i=1}^{m}{{w}_{i}}{{x}_{i}} \,|  $
 $ {\mathbf{x}} \geqslant X\mathbf{\lambda},\;{\mathbf{y}} \leqslant Y{\mathbf{\lambda },\;{\mathbf{\mathbf{e\lambda=1}, \lambda }} \geqslant {\mathbf{0}}} \Big\}$. Calculating maximum profit along with the optimal output and input quantities $\mathbf{y^{*}}$and $\mathbf{x^{*}}$ requires solving:

```math
\begin{aligned}
 & \underset{\mathbf{x,y,\lambda} }{\mathop{\max }}\,\quad \quad \quad \;\ \Pi\left(\mathbf{w},\mathbf{p}\right)=\mathbf{py^{*}-wx^{*}}  \\
 & \text{subject}\ \text{to} \\
 & \quad \quad \quad \quad \quad \ {{\mathbf{x}}}\ge X\mathbf{\lambda=x } \\
 & \quad \quad \quad \quad \quad  \; {{\mathbf{y}}}  \le Y\mathbf{\lambda =y} \\
 & \quad \quad \quad \quad \quad \; \mathbf{e\lambda=1} \\
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}.  
\end{aligned}
```

*Profit efficiency* defines as the difference between maximum profit and observed profit. Following the duality results introduced by *Chambers, Chung and Färe (1998)* it is possible to decompose it into technical and allocative efficiencies under variable returns to scale. Profit efficiency can be then decomposed into the directional distance fucntion and the residual difference corresponding to the *allocative profit efficiency*. Allocative efficiency defines then as the difference between maximum profit and profit at the technically efficient projection on the frontier. The approach relies on the directional vector to normalize these components, thereby ensuring that their values can be compared across DMUs.

In this example we compute the profit efficiency measure under variable returns to scale:
```@example revenue
using DataEnvelopmentAnalysis

X = [1 1; 1 1; 0.75 1.5; 0.5 2; 0.5 2; 2 2; 2.75 3.5; 1.375 1.75];

Y = [1 11; 5 3; 5 5; 2 9; 4 5; 4 2; 3 3; 4.5 3.5];

P = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

W = [2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1; 2 1];

deaprofit(X, Y, W, P, Gx = :Monetary, Gy = :Monetary)
```

### deaprofit Function Documentation

```@docs
deaprofit
```

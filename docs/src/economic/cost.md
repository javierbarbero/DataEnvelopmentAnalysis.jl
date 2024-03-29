```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Cost Models

## Cost Efficiency Model with Radial Technical Efficiency

Let us denote by $C\left(\mathbf{y},\mathbf{w}\right)$ the minimum cost of producing the output level $\mathbf{y}$ given the input price vector $\mathbf{w}$: $C\left(\mathbf{y},\mathbf{w}\right)=\min \left\{ \sum\limits_{i=1}^{m}{{{w}_{i}}{{x}_{i}}} | {\mathbf{x}} \geqslant X\mathbf{\lambda} {\mathbf{y}_{o}} \leqslant Y{\mathbf{\lambda },\;{\mathbf{\lambda }} \geqslant {\mathbf{0}}} \right\}$, which considers the input possibility set capable of producing $\mathbf{y}_{o}$. For the observed outputs levels we can calculate minimum cost and the associated optimal quantities of inputs $\mathbf{x^{*}}$ consistent with the production technology by solving the following program:

```math
\begin{aligned}
 & \underset{\mathbf{x} ,\mathbf{\lambda }}{\mathop{\min }}\,\quad \quad \quad \;\ C\left(\mathbf{y}_{},\mathbf{w}\right)=\mathbf{wx^{*}}  \\ 
 & \text{subject}\ \text{to} \\ 
 & \quad \quad \quad \quad \quad \ {{\mathbf{x}}}\ge X\mathbf{\lambda } \\ 
 & \quad \quad \quad \quad \quad  \;Y\mathbf{\lambda }\ \ge {{\mathbf{y}}_{o}}  \\ 
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. 
\end{aligned}
```

The measurement of cost efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

*Cost efficiency* defines as the ratio of minimum cost to observed cost: $CE=C\left(\mathbf{y},\mathbf{w}\right)/\mathbf{wx_{o}}$. Thanks to duality results presented by *Shephard (1953)* , and following *Farrell (1957)*, cost efficiency can be decomposed into the radially input oriented technical efficiency measure and the residual difference corresponding to *allocative cost efficiency*. Allocative  efficiency defines as the ratio between minimum cost to production cost at the technically efficient projection of the unit under evaluation.

In this example we compute the cost efficiency measure under variable returns to scale:
```@example cost
using DataEnvelopmentAnalysis

X = [5 3; 2 4; 4 2; 4 8; 7 9.0];

Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];

W = [2 1; 2 1; 2 1; 2 1; 2 1.0];

deacost(X, Y, W)
```

### deacost Function Documentation

```@docs
deacost
```


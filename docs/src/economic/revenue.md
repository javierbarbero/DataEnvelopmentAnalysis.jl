```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Revenue Models

## Revenue Efficiency Model with Radial Technical Efficiency

Let us denote by $R\left(\mathbf{x},\mathbf{p}\right)$ the maximum feasible revenue using inputs' levels $\mathbf{x}$ and given the outputs' prices $\mathbf{p}$: $R\left(\mathbf{x},\mathbf{p}\right)=\max \left\{ \sum\limits_{i=1}^{s}{{{p}_{i}}{{y}_{i}}} | {\mathbf{x}_{o}} \geqslant X\mathbf{\lambda},\;{\mathbf{y}} \leqslant Y{\mathbf{\lambda },\;{\mathbf{\lambda }} \geqslant {\mathbf{0}}} \right\}$; i.e.,  considering the output possibility set producible with $\mathbf{x}_{o}$. In this case, we calculate maximum revenue along with the optimal output quantities $\mathbf{y^{*}}$  by solving the following program:

```math
\begin{aligned}
 & \underset{\mathbf{y} ,\mathbf{\lambda }}{\mathop{\max }}\,\quad \quad \quad \;\ R\left(\mathbf{x}_{o},\mathbf{p}\right)=\mathbf{py^{*}}  \\ 
 & \text{subject}\ \text{to} \\ 
 & \quad \quad \quad \quad \quad \ {{\mathbf{x}_o}}\ge X\mathbf{\lambda } \\ 
 & \quad \quad \quad \quad \quad  \;Y\mathbf{\lambda }\ \ge {{\mathbf{y}}} \\ 
 & \quad \quad \quad \quad \quad \ \mathbf{\lambda }\ge \mathbf{0}. 
\end{aligned}
```

The measurement of revenue efficiency assuming variable returns to scale, **VRS**, adds the following condition:
```math
\sum\nolimits_{j=1}^{n}\lambda_j=1
```

*Revenue efficiency* defines as the ratio of observed revenue to maximum revenue: $RE=\mathbf{py_{o}}/R\left(\mathbf{x},\mathbf{p}\right) $. Duality results presented in *Shephard (1953)* from an output perspective allow us to decompose $RE$ into the output oriented technical efficiency measure and the residual difference corresponding to the *allocative revenue efficiency*. Allocative efficiency defines as the ratio of revenue at the technically efficient projection of the observation to maximum revenue.

In this example we compute the revnue efficiency measure under variable returns to scale:
```@example revenue
using DataEnvelopmentAnalysis

X = [5 3; 2 4; 4 2; 4 8; 7 9.0];

Y = [7 4; 10 8; 8 10; 5 4; 3 6.0];

P = [3 2; 3 2; 3 2; 3 2; 3 2.0];

dearevenue(X, Y, P)
```

### dearevenue Function Documentation

```@docs
dearevenue
```


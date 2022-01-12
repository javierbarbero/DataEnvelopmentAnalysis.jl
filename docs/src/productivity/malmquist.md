```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# The Malmquist index

## The Malmquist Productivity Index

The Malmquist index introduced by *Caves, Christensen and Diewert(1982)* measures the change in  productivity of the observation under evaluation by comparing its relative performance with respect to  reference  technologies corresponding to two different time periods.

Following *Fare, Grosskopf, Norris and Zhang (1994)* productivity change can be decomposed into efficiency change and technical change under the assumption of a constant returns to scale techncology.

In this example we compute the Malmquist productivity index:
```@example malmquist
using DataEnvelopmentAnalysis

X = Array{Float64,3}(undef, 5, 1, 2);
X[:, :, 1] = [2; 3; 5; 4; 4];
X[:, :, 2] = [1; 2; 4; 3; 4];

Y = Array{Float64,3}(undef, 5, 1, 2);
Y[:, :, 1] = [1; 4; 6; 3; 5];
Y[:, :, 2] = [1; 4; 6; 3; 3];

malmquist(X, Y)
```

### malmquist Function Documentation

```@docs
malmquist
```
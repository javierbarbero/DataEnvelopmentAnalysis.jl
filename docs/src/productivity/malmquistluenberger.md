```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# The Malmquist-Luenberger index

## The Malmquist-Luenberger Productivity Index

The Malmquist-Luenberger index (*Chung, Färe and Grosskopf, 1997*, and *Aparicio, Pastor and Zofío, 2013*) measures the change in  productivity of the observation under evaluation by comparing its relative performance with respect to  reference  technologies corresponding to two different time periods. Productivity change can be decomposed into efficiency change and technical change under the assumption of a constant returns to scale techncology.

In this example we compute the Malmquist-Luenberger productivity index for inputs `X`, good outputs `Y`, and bad outputs `B`, using directions `Gx`, `Gy`, and `Gb`. If not specified, default directions are `Gx = :Zeros`, `Gy = :Observed`, `Gb = :Observed`. If not specified, default directions are `Gx = :Zeros`, `Gy = :Observed`, `Gb = :Observed`.

```@example malmquistluenberger
using DataEnvelopmentAnalysis

X = Array{Float64,3}(undef, 5, 1, 2)
Y = Array{Float64,3}(undef, 5, 1, 2)
B = Array{Float64,3}(undef, 5, 1, 2)

X[:,:,1] = ones(5, 1);
X[:,:,2] = ones(5, 1);

Y[:,:,1] = [7; 5; 1; 3; 4];
Y[:,:,2] = [8; 5.5; 2; 2; 4];

B[:,:,1] = [2; 5; 3; 3; 2];
B[:,:,2] = [1; 3; 2; 4; 1];

mlprod = malmluen(X, Y, B)
```

### malmluen Function Documentation

```@docs
malmluen
```
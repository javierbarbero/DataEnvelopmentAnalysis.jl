```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Environmental Model

In this example we compute the *Aparicio, Pastor and Zofio (2023)* environmental model for inputs `X`, good outputs `Y`, and bad outputs `B`, using directions `Gx`, `Gy`, and `Gb`. If not specified, default directions are `Gx = :Zeros`, `Gy = :Observed`, `Gb = :Observed`.

```@example env
using DataEnvelopmentAnalysis

X = ones(5, 1);
Y = [7; 5; 1; 3; 4];
B = [2; 5; 3; 3; 2];

deaenv1 = deaenv(X, Y, B)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example env
efficiency(deaenv1)
```

The optimal peers, ``λ``, are returned with the `peers` function and are returned as a `DEAPeers` object:
```@example env
peers(deaenv1)
```

### deaenv Function Documentation

```@docs
deaenv
```

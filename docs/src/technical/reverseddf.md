```@meta
CurrentModule = DataEnvelopmentAnalysis
```

# Reverse Directional Distance Function

In this example, we compute the Reverse Directional Distance Function (Pastor et al., 2016) DEA model for the Enhanced Russell Graph associated efficiency measure under variable returns to scale:
```@example rddf
using DataEnvelopmentAnalysis

X = [2; 4; 8; 12; 6; 14; 14; 9.412];

Y = [1; 5; 8; 9; 3; 7; 9; 2.353];

dearddf(X, Y, :ERG, rts = :VRS)
```

Estimated efficiency scores are returned with the `efficiency` function:
```@example rddf
rddferg = dearddf(X, Y, :ERG, rts = :VRS);
nothing #hide
```

```@example rddf
efficiency(rddferg)
```

### dearddf Function Documentation

```@docs
dearddf
```

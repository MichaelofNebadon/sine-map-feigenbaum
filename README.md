# sine-map-feigenbaum
# Feigenbaum Universality in the Sine Map

A numerical verification of Feigenbaum universality for the sine map **x_{n+1} = r sin(π x_n)**, confirming that the universal constants δ and α converge to the Feigenbaum values to five to six significant figures using only NumPy.

-----

## Results

|Quantity                  |Computed (sine map)|Feigenbaum universal|Residual   |
|--------------------------|-------------------|--------------------|-----------|
|r* (accumulation point)   |0.8655792687       |0.8655792683 (ref.) |4.3 × 10⁻¹⁰|
|δ (parameter scaling, n=6)|4.668967           |4.669201609…        |2.3 × 10⁻⁴ |
|α (orbital scaling, n=6→7)|2.502943           |2.502907875…        |3.5 × 10⁻⁵ |

Both constants are confirmed to five to six significant figures. No stochastic averaging or heuristic fitting was used — results follow directly from the fixed-point structure of the map.

-----

## Quickstart

```bash
git clone https://github.com/MichaelofNebadon/sine-map-feigenbaum
cd sine-map-feigenbaum
python feigenbaum_sine_map.py
```

**Requirements:** Python 3, NumPy only.

-----

## Method

### Superstable cycle Newton solver

A 2ⁿ-cycle is *superstable* when it passes through the critical point x = 0.5, making the stability eigenvalue exactly zero. These cycles exist at unique parameter values rₙ located by Newton’s method applied to:

```
f^(2^n)(0.5; r) = 0.5
```

The total derivative dx/dr is tracked analytically along the orbit to form the Newton step:

```python
def superstable_r(n_doublings, r0):
    r = r0
    for _ in range(60):
        x, dxdr = 0.5, 0.0
        for _ in range(2 ** n_doublings):
            dxdr = sin(pi*x) + r*pi*cos(pi*x)*dxdr
            x = r * sin(pi * x)
        r -= (x - 0.5) / dxdr
    return r
```

### Feigenbaum δ

Successive ratios of parameter gaps converge to δ:

```
δₙ = (rₙ - rₙ₋₁) / (rₙ₊₁ - rₙ)  →  4.669201609…
```

### Feigenbaum α

The bifurcation splitting near x = 0.5 shrinks by α at each period-doubling. The second-nearest orbit point to x = 0.5 at the superstable parameter gives the splitting distance dₙ:

```
αₙ = dₙ / dₙ₊₁  →  2.502907875…
```

### Superstable cascade

|n|Cycle 2ⁿ|rₙ (superstable)|Gap rₙ - rₙ₋₁|
|-|--------|----------------|-------------|
|1|2       |0.7777337662    |—            |
|2|4       |0.8463821717    |0.0686484055 |
|3|8       |0.8614503509    |0.0150681792 |
|4|16      |0.8646941807    |0.0032438299 |
|5|32      |0.8653896734    |0.0006954927 |
|6|64      |0.8655386616    |0.0001489882 |
|7|128     |0.8655705719    |0.0000319103 |

### δ convergence

|n|δₙ      |Error from 4.6692016|
|-|--------|--------------------|
|2|4.555853|0.113349            |
|3|4.645182|0.024020            |
|4|4.664075|0.005127            |
|5|4.668106|0.001096            |
|6|4.668967|0.000235            |

### r* extrapolation

Using the geometric series tail `r* = rₙ + (rₙ - rₙ₋₁) / (δ - 1)`:

|From n|r* extrapolated|Error from target|
|------|---------------|-----------------|
|5     |0.8655792222   |4.6 × 10⁻⁸       |
|6     |0.8655792667   |1.6 × 10⁻⁹       |
|7     |0.8655792687   |4.3 × 10⁻¹⁰      |

### α convergence

|Transition n→n+1|αₙ      |Error from 2.5029079|
|----------------|--------|--------------------|
|1→2             |2.590697|0.087789            |
|2→3             |2.521384|0.018476            |
|3→4             |2.506683|0.003776            |
|4→5             |2.503703|0.000795            |
|5→6             |2.503075|0.000167            |
|6→7             |2.502943|0.000035            |

The error at each step shrinks by a factor of ~4.7 — itself equal to δ, as predicted by renormalization group theory.

-----

## Background

Feigenbaum’s 1978 discovery that δ and α are independent of the specific map function — depending only on the order of the critical point — is one of the deepest results in dynamical systems theory. The sine map and the logistic map `x → r x(1-x)` are different functions, but because both have a quadratic maximum, they share identical universal constants to arbitrary precision.

This repository verifies that result from first principles using only the superstable cycle structure of the map, with no external dependencies beyond NumPy.

-----

## Paper

A full write-up is available in [`Feigenbaum_Universality_Sine_Map.docx`](Feigenbaum_Universality_Sine_Map.docx).

-----

## License

MIT

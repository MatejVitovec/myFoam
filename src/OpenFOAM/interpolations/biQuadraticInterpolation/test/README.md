# biQuadraticInterpolation test

This test exercises coefficient generation, evaluation, all coordinate
transformations, and coefficient save/load round trips using `sin(x*y)` on
`[0.5, 1.5] x [0.5, 1.5]`. The mesh has two sections in each direction,
with 28 x-cells and 26 y-cells in total.

Build and run from this directory in an initialized OpenFOAM environment:

```sh
wmake
$FOAM_USER_APPBIN/testBiQuadraticInterpolation
```

The program writes one gnuplot data file and one coefficient dictionary for
each transform case in the current directory. The data columns are `x`, `y`,
exact value, computed value, reloaded value, interpolation error, and
save/load difference. It reports maximum and RMS interpolation error and
fails if the maximum save/load difference exceeds `1e-4`.

Plot the exact and interpolated surfaces with:

```sh
gnuplot -persist plotResults.gp
```

# biCubicInterpolation test

This test checks exact reproduction of a bilinear function for all coordinate
transformations and verifies coefficient save/load round trips. It also writes
sampled exact/interpolated values for `sin(x*y)` to gnuplot data files. The mesh
has two sections in each direction.

Build and run from this directory in an initialized OpenFOAM environment:

```sh
wmake
$FOAM_USER_APPBIN/testBiCubicInterpolation
```

The test passes when the maximum interpolation error is at most `1e-10` and
the maximum save/load difference is at most `1e-4`. It writes one `.dat` file
and one coefficient dictionary per transformation to the current directory.

Generate comparison and error images from that directory with:

```sh
gnuplot /path/to/biCubicInterpolation/test/plotResults.gp
```

This writes `test_<transformation>.png` images for all five transformation
cases.

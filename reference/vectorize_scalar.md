# Lift a scalar function to accept and return numeric vectors/arrays.

Scalar input gives a scalar `double` back; array input gives a `double`
array of the same shape. This is the direct analogue of Python's
`meijer.vectorize`, built on `vapply` rather than `base::Vectorize()` so
the scalar-in/scalar-out, array-in/same-shape-array-out contract is
enforced by the type signature, not merely likely.

## Usage

``` r
vectorize_scalar(scalar_fn)
```

## Arguments

  - scalar\_fn:
    
    A function taking a scalar (plus optional extra arguments) and
    returning a scalar.

## Value

A function with the same extra arguments as `scalar_fn`, taking `x` as a
scalar or array and returning a value of the same shape.

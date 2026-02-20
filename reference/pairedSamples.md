# View the relationship between miRNA and gene samples

This function allows to access the `pairedSamples` slot of a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object. The `MirnaExperiment` class is able to contain miRNA and gene
expression measurements deriving from the same individuals (paired
samples), or from different subjects (unpaired samples).

## Usage

``` r
pairedSamples(object)
```

## Arguments

- object:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

## Value

A `logical` value that is either `TRUE` for paired samples, or `FALSE`
for unpaired samples.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# check if an existing MirnaExperiment object derive from paired samples
pairedSamples(obj)
#> [1] TRUE
```

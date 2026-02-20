# Explore miRNA-target pairs

This function accesses the `targets` slot of a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object. After retrieving miRNA targets with the
[`getTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/getTargets.md)
function, the interactions between miRNAs and target genes are stored in
the `targets` slot and can be explored with this function.

## Usage

``` r
mirnaTargets(object)
```

## Arguments

- object:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

## Value

A `data.frame` object containing the interactions between miRNAs and
target genes, as retrieved with the
[`getTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/getTargets.md)
function.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# visualize targets
targets_df <- mirnaTargets(obj)
```

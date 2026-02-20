# Use custom miRNA-target interactions

Apart from the
[`getTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/getTargets.md)
function, users can use their own custom list of miRNA-target
interactions. This function allows to directly add custom interactions
to an existing
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object.

## Usage

``` r
setTargets(mirnaObj, tg)
```

## Arguments

- mirnaObj:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- tg:

  A `data.frame` object specifying miRNA-target interactions. Check the
  **details section** for the format

## Value

A
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object containing miRNA targets stored in the `targets` slot. Results
can be accessed with the
[`mirnaTargets()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaTargets.md)
function.

## Details

This function requires interactions to be passed as a `data.frame`
object with two column, one for miRNAs and one for genes. The column
names of this `data.frame` must include `MicroRNA` and `Gene.Symbol`.

## References

Ronchi, J., & Foti, M. (2026). MIRit: An integrative R framework for the
identification of impaired miRNA–mRNA regulatory networks in complex
diseases. Bioinformatics Advances, vbag042.
<https://doi.org/10.1093/bioadv/vbag042>

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# define targets
tg <- data.frame(MicroRNA = c("hsa-miR-1179", "hsa-miR-1179"),
                 Gene.Symbol = c("ACVR2A", "ABI2"))

# set targets
obj <- setTargets(obj, tg)

# print targets
mirnaTargets(obj)
#>       MicroRNA Gene.Symbol
#> 1 hsa-miR-1179      ACVR2A
#> 2 hsa-miR-1179        ABI2
```

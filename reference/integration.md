# Explore the results of the integration analysis between miRNAs and genes

After performing the integration analysis between miRNA and gene
expression values with the
[`mirnaIntegration()`](https://jacopo-ronchi.github.io/MIRit/reference/mirnaIntegration.md)
function, the results are stored in the `integration` slot of a
[`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
object and can be explored with this function.

## Usage

``` r
integration(object, param = FALSE)
```

## Arguments

- object:

  A
  [`MirnaExperiment`](https://jacopo-ronchi.github.io/MIRit/reference/MirnaExperiment-class.md)
  object containing miRNA and gene data

- param:

  Logical, whether to return the complete `list` object with the
  parameters used, or just the results stored in `data`. Default is
  FALSE

## Value

If `param` is FALSE, then this functions returns a `data.frame` object
containing the results of the integration analysis. Otherwise, a `list`
object including the parameters used for the analysis will be returned.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# load example MirnaExperiment object
obj <- loadExamples()

# visualize the results of correlation analysis
res <- integration(obj)
res
#>                             microRNA   Target microRNA.Direction
#> hsa.miR.1179.7          hsa-miR-1179    ITGA6      downregulated
#> hsa.miR.1179.10         hsa-miR-1179    RBMS2      downregulated
#> hsa.miR.1179.12         hsa-miR-1179   SPIRE1      downregulated
#> hsa.miR.1179.14         hsa-miR-1179     TAF5      downregulated
#> hsa.miR.1275.1          hsa-miR-1275    CCND2      downregulated
#> hsa.miR.1275.2          hsa-miR-1275     CD44      downregulated
#> hsa.miR.1275.7          hsa-miR-1275    ITGA3      downregulated
#> hsa.miR.1275.8          hsa-miR-1275    LLGL1      downregulated
#> hsa.miR.1275.9          hsa-miR-1275      MDK      downregulated
#> hsa.miR.1275.10         hsa-miR-1275     NAB2      downregulated
#> hsa.miR.1275.11         hsa-miR-1275     NUMB      downregulated
#> hsa.miR.1275.14         hsa-miR-1275     THRA      downregulated
#> hsa.miR.1275.16         hsa-miR-1275    TNNI1      downregulated
#> hsa.miR.138.1.3p    hsa-miR-138-1-3p   CCDC18      downregulated
#> hsa.miR.138.1.3p.1  hsa-miR-138-1-3p    CCND2      downregulated
#> hsa.miR.138.1.3p.8  hsa-miR-138-1-3p   LRP2BP      downregulated
#> hsa.miR.138.1.3p.11 hsa-miR-138-1-3p    TNNI1      downregulated
#> hsa.miR.138.1.3p.14 hsa-miR-138-1-3p   ZNF276      downregulated
#> hsa.miR.138.5p.1      hsa-miR-138-5p    ATG9B      downregulated
#> hsa.miR.138.5p.2      hsa-miR-138-5p     CD44      downregulated
#> hsa.miR.138.5p.3      hsa-miR-138-5p    CLCN2      downregulated
#> hsa.miR.138.5p.7      hsa-miR-138-5p     DTX4      downregulated
#> hsa.miR.138.5p.9      hsa-miR-138-5p     OGG1      downregulated
#> hsa.miR.138.5p.11     hsa-miR-138-5p    RBBP4      downregulated
#> hsa.miR.138.5p.14     hsa-miR-138-5p   VPS37A      downregulated
#> hsa.miR.139.5p.3      hsa-miR-139-5p    CCND2      downregulated
#> hsa.miR.139.5p.5      hsa-miR-139-5p     ECE1      downregulated
#> hsa.miR.139.5p.6      hsa-miR-139-5p      FN1      downregulated
#> hsa.miR.139.5p.8      hsa-miR-139-5p    ITGA3      downregulated
#> hsa.miR.139.5p.15     hsa-miR-139-5p   SPIRE1      downregulated
#> hsa.miR.139.5p.17     hsa-miR-139-5p   TSPAN3      downregulated
#> hsa.miR.144.3p.4      hsa-miR-144-3p   DCBLD2      downregulated
#> hsa.miR.144.3p.26     hsa-miR-144-3p   TSPAN3      downregulated
#> hsa.miR.144.5p.7      hsa-miR-144-5p     TAF5      downregulated
#> hsa.miR.146b.3p.1    hsa-miR-146b-3p     BCL2        upregulated
#> hsa.miR.146b.3p.9    hsa-miR-146b-3p    JMJD4        upregulated
#> hsa.miR.146b.3p.10   hsa-miR-146b-3p    LIMD1        upregulated
#> hsa.miR.146b.3p.11   hsa-miR-146b-3p   MED12L        upregulated
#> hsa.miR.146b.3p.13   hsa-miR-146b-3p     PAX8        upregulated
#> hsa.miR.146b.3p.15   hsa-miR-146b-3p   SORBS2        upregulated
#> hsa.miR.146b.3p.16   hsa-miR-146b-3p   SPTBN1        upregulated
#> hsa.miR.146b.3p.17   hsa-miR-146b-3p    TAF9B        upregulated
#> hsa.miR.146b.3p.22   hsa-miR-146b-3p  ZFP36L2        upregulated
#> hsa.miR.146b.3p.23   hsa-miR-146b-3p  ZMYND11        upregulated
#> hsa.miR.146b.5p      hsa-miR-146b-5p  ANKRD28        upregulated
#> hsa.miR.146b.5p.2    hsa-miR-146b-5p     CD24        upregulated
#> hsa.miR.146b.5p.3    hsa-miR-146b-5p     CUL5        upregulated
#> hsa.miR.146b.5p.4    hsa-miR-146b-5p   CYBRD1        upregulated
#> hsa.miR.146b.5p.5    hsa-miR-146b-5p     DIO2        upregulated
#> hsa.miR.146b.5p.6    hsa-miR-146b-5p     EGR1        upregulated
#> hsa.miR.146b.5p.7    hsa-miR-146b-5p    FBXL3        upregulated
#> hsa.miR.146b.5p.8    hsa-miR-146b-5p     FHL1        upregulated
#> hsa.miR.146b.5p.10   hsa-miR-146b-5p    HYOU1        upregulated
#> hsa.miR.146b.5p.15   hsa-miR-146b-5p     PAX8        upregulated
#> hsa.miR.146b.5p.16   hsa-miR-146b-5p    SMAD2        upregulated
#> hsa.miR.146b.5p.17   hsa-miR-146b-5p   SPTBN1        upregulated
#> hsa.miR.146b.5p.18   hsa-miR-146b-5p    TAF9B        upregulated
#> hsa.miR.182.5p        hsa-miR-182-5p     BCL2        upregulated
#> hsa.miR.182.5p.2      hsa-miR-182-5p   CITED2        upregulated
#> hsa.miR.182.5p.3      hsa-miR-182-5p     CUL5        upregulated
#> hsa.miR.182.5p.4      hsa-miR-182-5p     FHL1        upregulated
#> hsa.miR.182.5p.8      hsa-miR-182-5p     SDC2        upregulated
#> hsa.miR.182.5p.10     hsa-miR-182-5p   SPTBN1        upregulated
#> hsa.miR.182.5p.12     hsa-miR-182-5p     TOB1        upregulated
#> hsa.miR.183.5p.6      hsa-miR-183-5p      EZR        upregulated
#> hsa.miR.183.5p.7      hsa-miR-183-5p   HNRNPM        upregulated
#> hsa.miR.183.5p.11     hsa-miR-183-5p   MED12L        upregulated
#> hsa.miR.204.5p.1      hsa-miR-204-5p    CCND2      downregulated
#> hsa.miR.204.5p.2      hsa-miR-204-5p     CD44      downregulated
#> hsa.miR.204.5p.3      hsa-miR-204-5p     DTX4      downregulated
#> hsa.miR.204.5p.8      hsa-miR-204-5p HLA-DRB5      downregulated
#> hsa.miR.204.5p.10     hsa-miR-204-5p  RHOBTB2      downregulated
#> hsa.miR.204.5p.11     hsa-miR-204-5p     SLA2      downregulated
#> hsa.miR.204.5p.13     hsa-miR-204-5p    STAT6      downregulated
#> hsa.miR.204.5p.14     hsa-miR-204-5p     TAF5      downregulated
#> hsa.miR.21.3p          hsa-miR-21-3p     CANX        upregulated
#> hsa.miR.21.3p.1        hsa-miR-21-3p     CUL5        upregulated
#> hsa.miR.21.3p.2        hsa-miR-21-3p     DIO2        upregulated
#> hsa.miR.21.3p.3        hsa-miR-21-3p     GJA1        upregulated
#> hsa.miR.21.3p.7        hsa-miR-21-3p     SDC2        upregulated
#> hsa.miR.21.3p.8        hsa-miR-21-3p  SMARCA2        upregulated
#> hsa.miR.21.5p          hsa-miR-21-5p    AIF1L        upregulated
#> hsa.miR.21.5p.1        hsa-miR-21-5p  ANKRD28        upregulated
#> hsa.miR.21.5p.3        hsa-miR-21-5p    APAF1        upregulated
#> hsa.miR.21.5p.4        hsa-miR-21-5p ARHGAP24        upregulated
#> hsa.miR.21.5p.5        hsa-miR-21-5p     BCL2        upregulated
#> hsa.miR.21.5p.6        hsa-miR-21-5p   CYBRD1        upregulated
#> hsa.miR.21.5p.8        hsa-miR-21-5p    MATN2        upregulated
#> hsa.miR.21.5p.10       hsa-miR-21-5p   MED12L        upregulated
#> hsa.miR.21.5p.15       hsa-miR-21-5p     TOB1        upregulated
#> hsa.miR.2110.1          hsa-miR-2110    AIF1L        upregulated
#> hsa.miR.2110.3          hsa-miR-2110   CYBRD1        upregulated
#> hsa.miR.2110.9          hsa-miR-2110   MED12L        upregulated
#> hsa.miR.2110.12         hsa-miR-2110    SMAD2        upregulated
#> hsa.miR.221.3p        hsa-miR-221-3p  ANKRD28        upregulated
#> hsa.miR.221.3p.2      hsa-miR-221-3p    APAF1        upregulated
#> hsa.miR.221.3p.3      hsa-miR-221-3p     BCL2        upregulated
#> hsa.miR.221.3p.6      hsa-miR-221-3p     CUL5        upregulated
#> hsa.miR.221.3p.17     hsa-miR-221-3p  SMARCA2        upregulated
#> hsa.miR.221.3p.18     hsa-miR-221-3p   SPTBN1        upregulated
#> hsa.miR.221.3p.21     hsa-miR-221-3p     URI1        upregulated
#> hsa.miR.221.5p        hsa-miR-221-5p  ANKRD28        upregulated
#> hsa.miR.221.5p.1      hsa-miR-221-5p    APAF1        upregulated
#> hsa.miR.221.5p.6      hsa-miR-221-5p      IYD        upregulated
#> hsa.miR.221.5p.7      hsa-miR-221-5p    LIMD1        upregulated
#> hsa.miR.221.5p.10     hsa-miR-221-5p     URI1        upregulated
#> hsa.miR.222.3p        hsa-miR-222-3p  ANKRD28        upregulated
#> hsa.miR.222.3p.2      hsa-miR-222-3p    APAF1        upregulated
#> hsa.miR.222.3p.4      hsa-miR-222-3p     BCL2        upregulated
#> hsa.miR.222.3p.5      hsa-miR-222-3p     CANX        upregulated
#> hsa.miR.222.3p.9      hsa-miR-222-3p      EZR        upregulated
#> hsa.miR.222.3p.11     hsa-miR-222-3p    GPBP1        upregulated
#> hsa.miR.222.3p.18     hsa-miR-222-3p  SMARCA2        upregulated
#> hsa.miR.222.3p.22     hsa-miR-222-3p     URI1        upregulated
#> hsa.miR.31.3p          hsa-miR-31-3p     CANX        upregulated
#> hsa.miR.31.3p.1        hsa-miR-31-3p   COL9A3        upregulated
#> hsa.miR.31.3p.2        hsa-miR-31-3p     CUL5        upregulated
#> hsa.miR.31.3p.3        hsa-miR-31-3p    LIMD1        upregulated
#> hsa.miR.31.3p.4        hsa-miR-31-3p   MED12L        upregulated
#> hsa.miR.31.3p.6        hsa-miR-31-3p  SLC26A7        upregulated
#> hsa.miR.31.3p.7        hsa-miR-31-3p    SMAD2        upregulated
#> hsa.miR.31.5p.2        hsa-miR-31-5p     CUL5        upregulated
#> hsa.miR.31.5p.5        hsa-miR-31-5p   GIGYF2        upregulated
#> hsa.miR.31.5p.6        hsa-miR-31-5p      JUN        upregulated
#> hsa.miR.31.5p.7        hsa-miR-31-5p     MICA        upregulated
#> hsa.miR.31.5p.10       hsa-miR-31-5p   SPTBN1        upregulated
#> hsa.miR.34a.5p        hsa-miR-34a-5p    APAF1        upregulated
#> hsa.miR.34a.5p.1      hsa-miR-34a-5p     BCL2        upregulated
#> hsa.miR.34a.5p.3      hsa-miR-34a-5p     CD24        upregulated
#> hsa.miR.34a.5p.7      hsa-miR-34a-5p   GIGYF2        upregulated
#> hsa.miR.34a.5p.9      hsa-miR-34a-5p    LIMD1        upregulated
#> hsa.miR.34a.5p.13     hsa-miR-34a-5p   SPTBN1        upregulated
#> hsa.miR.34a.5p.14     hsa-miR-34a-5p    SYNE2        upregulated
#> hsa.miR.34a.5p.18     hsa-miR-34a-5p  ZMYND11        upregulated
#> hsa.miR.3613.5p.3    hsa-miR-3613-5p   MED12L        upregulated
#> hsa.miR.3613.5p.4    hsa-miR-3613-5p   PAPSS2        upregulated
#> hsa.miR.3613.5p.5    hsa-miR-3613-5p  SLC26A4        upregulated
#> hsa.miR.3613.5p.6    hsa-miR-3613-5p  SLC26A7        upregulated
#> hsa.miR.3613.5p.7    hsa-miR-3613-5p    SYNE2        upregulated
#> hsa.miR.375.1            hsa-miR-375     BCL2        upregulated
#> hsa.miR.375.3            hsa-miR-375      DEK        upregulated
#> hsa.miR.451a.2          hsa-miR-451a    CCND2      downregulated
#> hsa.miR.451a.3          hsa-miR-451a   DCBLD2      downregulated
#> hsa.miR.451a.7          hsa-miR-451a   TNRC6C      downregulated
#> hsa.miR.486.3p.3      hsa-miR-486-3p    CCND2      downregulated
#> hsa.miR.486.3p.5      hsa-miR-486-3p    EFHD2      downregulated
#> hsa.miR.486.3p.9      hsa-miR-486-3p    LLGL1      downregulated
#> hsa.miR.486.5p        hsa-miR-486-5p   DCBLD2      downregulated
#> hsa.miR.486.5p.5      hsa-miR-486-5p    RBMS2      downregulated
#> hsa.miR.504.5p.6      hsa-miR-504-5p    LLGL1      downregulated
#> hsa.miR.504.5p.7      hsa-miR-504-5p    LUC7L      downregulated
#> hsa.miR.504.5p.9      hsa-miR-504-5p   MTHFD2      downregulated
#> hsa.miR.504.5p.10     hsa-miR-504-5p     NAB2      downregulated
#> hsa.miR.504.5p.11     hsa-miR-504-5p    TNNI1      downregulated
#> hsa.miR.551b.3p      hsa-miR-551b-3p    AIF1L        upregulated
#> hsa.miR.551b.3p.3    hsa-miR-551b-3p   GIGYF2        upregulated
#> hsa.miR.551b.3p.4    hsa-miR-551b-3p    GPBP1        upregulated
#> hsa.miR.551b.3p.8    hsa-miR-551b-3p      IYD        upregulated
#> hsa.miR.551b.3p.10   hsa-miR-551b-3p    NUPR1        upregulated
#> hsa.miR.551b.3p.12   hsa-miR-551b-3p     RYR2        upregulated
#> hsa.miR.551b.3p.13   hsa-miR-551b-3p  SMARCA2        upregulated
#> hsa.miR.551b.3p.14   hsa-miR-551b-3p    SMPD4        upregulated
#> hsa.miR.577              hsa-miR-577      AGK      downregulated
#> hsa.miR.577.4            hsa-miR-577     CD44      downregulated
#> hsa.miR.577.8            hsa-miR-577     LMO3      downregulated
#> hsa.miR.577.9            hsa-miR-577   LRP2BP      downregulated
#> hsa.miR.577.10           hsa-miR-577   MTHFD2      downregulated
#> hsa.miR.577.17           hsa-miR-577     WSB2      downregulated
#> hsa.miR.652.3p.4      hsa-miR-652-3p    CCND2      downregulated
#> hsa.miR.652.3p.6      hsa-miR-652-3p    EFHD2      downregulated
#> hsa.miR.652.3p.11     hsa-miR-652-3p    LLGL1      downregulated
#> hsa.miR.652.3p.12     hsa-miR-652-3p     LMO3      downregulated
#> hsa.miR.652.3p.13     hsa-miR-652-3p    MGST1      downregulated
#> hsa.miR.652.3p.18     hsa-miR-652-3p   VPS37A      downregulated
#> hsa.miR.652.3p.20     hsa-miR-652-3p  ZFP36L1      downregulated
#> hsa.miR.653.5p.1      hsa-miR-653-5p     BCL2        upregulated
#> hsa.miR.653.5p.3      hsa-miR-653-5p   CITED2        upregulated
#> hsa.miR.653.5p.6      hsa-miR-653-5p     DIO2        upregulated
#> hsa.miR.653.5p.10     hsa-miR-653-5p    HYOU1        upregulated
#> hsa.miR.653.5p.12     hsa-miR-653-5p      IYD        upregulated
#> hsa.miR.653.5p.15     hsa-miR-653-5p   NECAB1        upregulated
#> hsa.miR.653.5p.21     hsa-miR-653-5p    SMAD2        upregulated
#> hsa.miR.653.5p.22     hsa-miR-653-5p     URI1        upregulated
#> hsa.miR.6842.3p.4    hsa-miR-6842-3p    LIMD1        upregulated
#> hsa.miR.6842.3p.7    hsa-miR-6842-3p  SLC26A7        upregulated
#> hsa.miR.7.5p.6          hsa-miR-7-5p    CCND2      downregulated
#> hsa.miR.873.3p        hsa-miR-873-3p    CCND2      downregulated
#> hsa.miR.873.3p.1      hsa-miR-873-3p     CD44      downregulated
#> hsa.miR.873.3p.7      hsa-miR-873-3p    RBBP4      downregulated
#> hsa.miR.873.3p.10     hsa-miR-873-3p    SYTL4      downregulated
#> hsa.miR.873.3p.12     hsa-miR-873-3p   VPS37A      downregulated
#> hsa.miR.873.5p.1      hsa-miR-873-5p     DDR1      downregulated
#> hsa.miR.873.5p.7      hsa-miR-873-5p   NPEPPS      downregulated
#> hsa.miR.873.5p.11     hsa-miR-873-5p    RBMS2      downregulated
#> hsa.miR.874.3p.5      hsa-miR-874-3p    CCND2      downregulated
#> hsa.miR.874.3p.6      hsa-miR-874-3p     ECE1      downregulated
#> hsa.miR.874.3p.7      hsa-miR-874-3p     FZD4      downregulated
#> hsa.miR.874.3p.11     hsa-miR-874-3p  S100A13      downregulated
#> hsa.miR.874.3p.15     hsa-miR-874-3p     THRA      downregulated
#> hsa.miR.874.3p.16     hsa-miR-874-3p    TNNI1      downregulated
#> hsa.miR.874.3p.17     hsa-miR-874-3p   TNRC6C      downregulated
#> hsa.miR.9.5p.2          hsa-miR-9-5p  C1orf56      downregulated
#> hsa.miR.9.5p.3          hsa-miR-9-5p   DCBLD2      downregulated
#> hsa.miR.9.5p.9          hsa-miR-9-5p   MTHFD2      downregulated
#> hsa.miR.9.5p.12         hsa-miR-9-5p     SIX5      downregulated
#> hsa.miR.96.5p.1        hsa-miR-96-5p ARHGAP24        upregulated
#> hsa.miR.96.5p.2        hsa-miR-96-5p     BCL2        upregulated
#> hsa.miR.96.5p.3        hsa-miR-96-5p     CAV1        upregulated
#> hsa.miR.96.5p.5        hsa-miR-96-5p      DEK        upregulated
#> hsa.miR.96.5p.8        hsa-miR-96-5p      EZR        upregulated
#> hsa.miR.96.5p.9        hsa-miR-96-5p     FHL1        upregulated
#> hsa.miR.96.5p.15       hsa-miR-96-5p   MED12L        upregulated
#> hsa.miR.96.5p.18       hsa-miR-96-5p     SDC2        upregulated
#> hsa.miR.96.5p.19       hsa-miR-96-5p  SLC26A4        upregulated
#> hsa.miR.96.5p.20       hsa-miR-96-5p    SYNE2        upregulated
#>                     Corr.Coefficient Corr.P.Value Corr.Adjusted.P.Val
#> hsa.miR.1179.7            -0.8617647 8.902592e-06        6.579303e-04
#> hsa.miR.1179.10           -0.7764706 2.021754e-04        3.484055e-03
#> hsa.miR.1179.12           -0.7117647 9.920216e-04        7.248522e-03
#> hsa.miR.1179.14           -0.5970588 7.304866e-03        2.668766e-02
#> hsa.miR.1275.1            -0.7794118 1.858401e-04        3.484055e-03
#> hsa.miR.1275.2            -0.8382353 2.505171e-05        1.342354e-03
#> hsa.miR.1275.7            -0.5647059 1.133160e-02        3.755784e-02
#> hsa.miR.1275.8            -0.6676471 2.355991e-03        1.165309e-02
#> hsa.miR.1275.9            -0.6823529 1.793985e-03        9.944246e-03
#> hsa.miR.1275.10           -0.7264706 7.180913e-04        5.996529e-03
#> hsa.miR.1275.11           -0.6705882 2.233604e-03        1.130872e-02
#> hsa.miR.1275.14           -0.5735294 1.009552e-02        3.490010e-02
#> hsa.miR.1275.16           -0.6852941 1.695781e-03        9.649443e-03
#> hsa.miR.138.1.3p          -0.6235294 4.927825e-03        2.044253e-02
#> hsa.miR.138.1.3p.1        -0.8941176 1.505655e-06        2.420340e-04
#> hsa.miR.138.1.3p.8        -0.5529412 1.315725e-02        4.147114e-02
#> hsa.miR.138.1.3p.11       -0.7764706 2.021754e-04        3.484055e-03
#> hsa.miR.138.1.3p.14       -0.5676471 1.090730e-02        3.691259e-02
#> hsa.miR.138.5p.1          -0.5647059 1.133160e-02        3.755784e-02
#> hsa.miR.138.5p.2          -0.7176471 8.737675e-04        6.688482e-03
#> hsa.miR.138.5p.3          -0.7029412 1.193476e-03        7.993805e-03
#> hsa.miR.138.5p.7          -0.7441176 4.741414e-04        4.917305e-03
#> hsa.miR.138.5p.9          -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.138.5p.11         -0.5970588 7.304866e-03        2.668766e-02
#> hsa.miR.138.5p.14         -0.6735294 2.116370e-03        1.080020e-02
#> hsa.miR.139.5p.3          -0.5588235 1.221826e-02        3.967848e-02
#> hsa.miR.139.5p.5          -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.139.5p.6          -0.6970588 1.345295e-03        8.317543e-03
#> hsa.miR.139.5p.8          -0.6794118 1.896716e-03        1.042383e-02
#> hsa.miR.139.5p.15         -0.5705882 1.049535e-02        3.608828e-02
#> hsa.miR.139.5p.17         -0.5411765 1.520011e-02        4.654129e-02
#> hsa.miR.144.3p.4          -0.6852941 1.695781e-03        9.649443e-03
#> hsa.miR.144.3p.26         -0.6147059 5.639737e-03        2.184549e-02
#> hsa.miR.144.5p.7          -0.7058824 1.122950e-03        7.934693e-03
#> hsa.miR.146b.3p.1         -0.8294118 3.544227e-05        1.498483e-03
#> hsa.miR.146b.3p.9         -0.6823529 1.793985e-03        9.944246e-03
#> hsa.miR.146b.3p.10        -0.6911765 1.512359e-03        8.840428e-03
#> hsa.miR.146b.3p.11        -0.6470588 3.371039e-03        1.548270e-02
#> hsa.miR.146b.3p.13        -0.6764706 2.004126e-03        1.047685e-02
#> hsa.miR.146b.3p.15        -0.8264706 3.961774e-05        1.498483e-03
#> hsa.miR.146b.3p.16        -0.8500000 1.526130e-05        8.920924e-04
#> hsa.miR.146b.3p.17        -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.146b.3p.22        -0.5911765 7.937790e-03        2.835555e-02
#> hsa.miR.146b.3p.23        -0.7264706 7.180913e-04        5.996529e-03
#> hsa.miR.146b.5p           -0.8235294 4.419570e-05        1.578769e-03
#> hsa.miR.146b.5p.2         -0.7117647 9.920216e-04        7.248522e-03
#> hsa.miR.146b.5p.3         -0.6852941 1.695781e-03        9.649443e-03
#> hsa.miR.146b.5p.4         -0.8029412 9.025840e-05        2.418173e-03
#> hsa.miR.146b.5p.5         -0.7029412 1.193476e-03        7.993805e-03
#> hsa.miR.146b.5p.6         -0.5382353 1.574661e-02        4.775978e-02
#> hsa.miR.146b.5p.7         -0.5529412 1.315725e-02        4.147114e-02
#> hsa.miR.146b.5p.8         -0.7176471 8.737675e-04        6.688482e-03
#> hsa.miR.146b.5p.10        -0.8264706 3.961774e-05        1.498483e-03
#> hsa.miR.146b.5p.15        -0.7117647 9.920216e-04        7.248522e-03
#> hsa.miR.146b.5p.16        -0.5352941 1.630793e-02        4.877209e-02
#> hsa.miR.146b.5p.17        -0.8352941 2.818566e-05        1.394106e-03
#> hsa.miR.146b.5p.18        -0.6176471 5.393992e-03        2.114840e-02
#> hsa.miR.182.5p            -0.8323529 3.164062e-05        1.453209e-03
#> hsa.miR.182.5p.2          -0.6676471 2.355991e-03        1.165309e-02
#> hsa.miR.182.5p.3          -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.182.5p.4          -0.6411765 3.717210e-03        1.671445e-02
#> hsa.miR.182.5p.8          -0.5617647 1.176851e-02        3.860791e-02
#> hsa.miR.182.5p.10         -0.5352941 1.630793e-02        4.877209e-02
#> hsa.miR.182.5p.12         -0.7205882 8.190956e-04        6.502203e-03
#> hsa.miR.183.5p.6          -0.5882353 8.269861e-03        2.921715e-02
#> hsa.miR.183.5p.7          -0.6294118 4.494045e-03        1.939376e-02
#> hsa.miR.183.5p.11         -0.6823529 1.793985e-03        9.944246e-03
#> hsa.miR.204.5p.1          -0.8852941 2.573574e-06        3.266236e-04
#> hsa.miR.204.5p.2          -0.8088235 7.423780e-05        2.273091e-03
#> hsa.miR.204.5p.3          -0.7588235 3.270151e-04        4.043668e-03
#> hsa.miR.204.5p.8          -0.5941176 7.616200e-03        2.751245e-02
#> hsa.miR.204.5p.10         -0.9029412 8.388622e-07        2.420340e-04
#> hsa.miR.204.5p.11         -0.7205882 8.190956e-04        6.502203e-03
#> hsa.miR.204.5p.13         -0.8205882 4.920635e-05        1.665247e-03
#> hsa.miR.204.5p.14         -0.7205882 8.190956e-04        6.502203e-03
#> hsa.miR.21.3p             -0.6617647 2.616877e-03        1.265152e-02
#> hsa.miR.21.3p.1           -0.7558824 3.529453e-04        4.281959e-03
#> hsa.miR.21.3p.2           -0.7647059 2.798381e-04        3.911651e-03
#> hsa.miR.21.3p.3           -0.5911765 7.937790e-03        2.835555e-02
#> hsa.miR.21.3p.7           -0.8617647 8.902592e-06        6.579303e-04
#> hsa.miR.21.3p.8           -0.6764706 2.004126e-03        1.047685e-02
#> hsa.miR.21.5p             -0.7882353 1.432328e-04        3.076035e-03
#> hsa.miR.21.5p.1           -0.7970588 1.090490e-04        2.696867e-03
#> hsa.miR.21.5p.3           -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.21.5p.4           -0.7882353 1.432328e-04        3.076035e-03
#> hsa.miR.21.5p.5           -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.21.5p.6           -0.6764706 2.004126e-03        1.047685e-02
#> hsa.miR.21.5p.8           -0.9352941 5.360084e-08        3.446534e-05
#> hsa.miR.21.5p.10          -0.7323529 6.274791e-04        5.682663e-03
#> hsa.miR.21.5p.15          -0.7617647 3.026714e-04        4.043668e-03
#> hsa.miR.2110.1            -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.2110.3            -0.6941176 1.426853e-03        8.574453e-03
#> hsa.miR.2110.9            -0.7823529 1.706104e-04        3.324317e-03
#> hsa.miR.2110.12           -0.6647059 2.483693e-03        1.209860e-02
#> hsa.miR.221.3p            -0.6411765 3.717210e-03        1.671445e-02
#> hsa.miR.221.3p.2          -0.6970588 1.345295e-03        8.317543e-03
#> hsa.miR.221.3p.3          -0.7441176 4.741414e-04        4.917305e-03
#> hsa.miR.221.3p.6          -0.7617647 3.026714e-04        4.043668e-03
#> hsa.miR.221.3p.17         -0.7382353 5.464254e-04        5.323508e-03
#> hsa.miR.221.3p.18         -0.5500000 1.364699e-02        4.239136e-02
#> hsa.miR.221.3p.21         -0.7382353 5.464254e-04        5.323508e-03
#> hsa.miR.221.5p            -0.5970588 7.304866e-03        2.668766e-02
#> hsa.miR.221.5p.1          -0.6911765 1.512359e-03        8.840428e-03
#> hsa.miR.221.5p.6          -0.7264706 7.180913e-04        5.996529e-03
#> hsa.miR.221.5p.7          -0.7941176 1.195955e-04        2.848144e-03
#> hsa.miR.221.5p.10         -0.6970588 1.345295e-03        8.317543e-03
#> hsa.miR.222.3p            -0.6323529 4.288792e-03        1.863306e-02
#> hsa.miR.222.3p.2          -0.7000000 1.267547e-03        8.316662e-03
#> hsa.miR.222.3p.4          -0.7058824 1.122950e-03        7.934693e-03
#> hsa.miR.222.3p.5          -0.7411765 5.092348e-04        5.197428e-03
#> hsa.miR.222.3p.9          -0.6970588 1.345295e-03        8.317543e-03
#> hsa.miR.222.3p.11         -0.5500000 1.364699e-02        4.239136e-02
#> hsa.miR.222.3p.18         -0.7323529 6.274791e-04        5.682663e-03
#> hsa.miR.222.3p.22         -0.7176471 8.737675e-04        6.688482e-03
#> hsa.miR.31.3p             -0.5382353 1.574661e-02        4.775978e-02
#> hsa.miR.31.3p.1           -0.5676471 1.090730e-02        3.691259e-02
#> hsa.miR.31.3p.2           -0.6235294 4.927825e-03        2.044253e-02
#> hsa.miR.31.3p.3           -0.6558824 2.900368e-03        1.381435e-02
#> hsa.miR.31.3p.4           -0.6500000 3.207853e-03        1.505584e-02
#> hsa.miR.31.3p.6           -0.6941176 1.426853e-03        8.574453e-03
#> hsa.miR.31.3p.7           -0.5558824 1.268109e-02        4.056686e-02
#> hsa.miR.31.5p.2           -0.6088235 6.157602e-03        2.342803e-02
#> hsa.miR.31.5p.5           -0.6029412 6.712063e-03        2.523893e-02
#> hsa.miR.31.5p.6           -0.6470588 3.371039e-03        1.548270e-02
#> hsa.miR.31.5p.7           -0.6764706 2.004126e-03        1.047685e-02
#> hsa.miR.31.5p.10          -0.7294118 6.715393e-04        5.915066e-03
#> hsa.miR.34a.5p            -0.7970588 1.090490e-04        2.696867e-03
#> hsa.miR.34a.5p.1          -0.8588235 1.023220e-05        6.579303e-04
#> hsa.miR.34a.5p.3          -0.5529412 1.315725e-02        4.147114e-02
#> hsa.miR.34a.5p.7          -0.5941176 7.616200e-03        2.751245e-02
#> hsa.miR.34a.5p.9          -0.6235294 4.927825e-03        2.044253e-02
#> hsa.miR.34a.5p.13         -0.6911765 1.512359e-03        8.840428e-03
#> hsa.miR.34a.5p.14         -0.6264706 4.706989e-03        1.991180e-02
#> hsa.miR.34a.5p.18         -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.3613.5p.3         -0.8823529 3.047810e-06        3.266236e-04
#> hsa.miR.3613.5p.4         -0.7823529 1.706104e-04        3.324317e-03
#> hsa.miR.3613.5p.5         -0.6970588 1.345295e-03        8.317543e-03
#> hsa.miR.3613.5p.6         -0.6941176 1.426853e-03        8.574453e-03
#> hsa.miR.3613.5p.7         -0.5647059 1.133160e-02        3.755784e-02
#> hsa.miR.375.1             -0.5617647 1.176851e-02        3.860791e-02
#> hsa.miR.375.3             -0.7323529 6.274791e-04        5.682663e-03
#> hsa.miR.451a.2            -0.8088235 7.423780e-05        2.273091e-03
#> hsa.miR.451a.3            -0.6676471 2.355991e-03        1.165309e-02
#> hsa.miR.451a.7            -0.6117647 5.894203e-03        2.269444e-02
#> hsa.miR.486.3p.3          -0.5588235 1.221826e-02        3.967848e-02
#> hsa.miR.486.3p.5          -0.6735294 2.116370e-03        1.080020e-02
#> hsa.miR.486.3p.9          -0.7088235 1.055841e-03        7.628152e-03
#> hsa.miR.486.5p            -0.6323529 4.288792e-03        1.863306e-02
#> hsa.miR.486.5p.5          -0.5676471 1.090730e-02        3.691259e-02
#> hsa.miR.504.5p.6          -0.7588235 3.270151e-04        4.043668e-03
#> hsa.miR.504.5p.7          -0.5441176 1.466817e-02        4.534439e-02
#> hsa.miR.504.5p.9          -0.5852941 8.612643e-03        3.026191e-02
#> hsa.miR.504.5p.10         -0.6264706 4.706989e-03        1.991180e-02
#> hsa.miR.504.5p.11         -0.5970588 7.304866e-03        2.668766e-02
#> hsa.miR.551b.3p           -0.7382353 5.464254e-04        5.323508e-03
#> hsa.miR.551b.3p.3         -0.5735294 1.009552e-02        3.490010e-02
#> hsa.miR.551b.3p.4         -0.8588235 1.023220e-05        6.579303e-04
#> hsa.miR.551b.3p.8         -0.7000000 1.267547e-03        8.316662e-03
#> hsa.miR.551b.3p.10        -0.7235294 7.672406e-04        6.324817e-03
#> hsa.miR.551b.3p.12        -0.7881698 1.435164e-04        3.076035e-03
#> hsa.miR.551b.3p.13        -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.551b.3p.14        -0.7676471 2.584418e-04        3.692846e-03
#> hsa.miR.577               -0.6352941 4.091032e-03        1.814161e-02
#> hsa.miR.577.4             -0.6470588 3.371039e-03        1.548270e-02
#> hsa.miR.577.8             -0.6764706 2.004126e-03        1.047685e-02
#> hsa.miR.577.9             -0.7029412 1.193476e-03        7.993805e-03
#> hsa.miR.577.10            -0.5500000 1.364699e-02        4.239136e-02
#> hsa.miR.577.17            -0.6558824 2.900368e-03        1.381435e-02
#> hsa.miR.652.3p.4          -0.8058824 8.192333e-05        2.394396e-03
#> hsa.miR.652.3p.6          -0.6500000 3.207853e-03        1.505584e-02
#> hsa.miR.652.3p.11         -0.7264706 7.180913e-04        5.996529e-03
#> hsa.miR.652.3p.12         -0.7117647 9.920216e-04        7.248522e-03
#> hsa.miR.652.3p.13         -0.7441176 4.741414e-04        4.917305e-03
#> hsa.miR.652.3p.18         -0.6147059 5.639737e-03        2.184549e-02
#> hsa.miR.652.3p.20         -0.6411765 3.717210e-03        1.671445e-02
#> hsa.miR.653.5p.1          -0.7470588 4.410535e-04        4.889610e-03
#> hsa.miR.653.5p.3          -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.653.5p.6          -0.6323529 4.288792e-03        1.863306e-02
#> hsa.miR.653.5p.10         -0.7588235 3.270151e-04        4.043668e-03
#> hsa.miR.653.5p.12         -0.6647059 2.483693e-03        1.209860e-02
#> hsa.miR.653.5p.15         -0.6764706 2.004126e-03        1.047685e-02
#> hsa.miR.653.5p.21         -0.5352941 1.630793e-02        4.877209e-02
#> hsa.miR.653.5p.22         -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.6842.3p.4         -0.7441176 4.741414e-04        4.917305e-03
#> hsa.miR.6842.3p.7         -0.5411765 1.520011e-02        4.654129e-02
#> hsa.miR.7.5p.6            -0.6735294 2.116370e-03        1.080020e-02
#> hsa.miR.873.3p            -0.8970588 1.246188e-06        2.420340e-04
#> hsa.miR.873.3p.1          -0.7323529 6.274791e-04        5.682663e-03
#> hsa.miR.873.3p.7          -0.7852941 1.564270e-04        3.244598e-03
#> hsa.miR.873.3p.10         -0.7352941 5.858078e-04        5.622006e-03
#> hsa.miR.873.3p.12         -0.6176471 5.393992e-03        2.114840e-02
#> hsa.miR.873.5p.1          -0.7294118 6.715393e-04        5.915066e-03
#> hsa.miR.873.5p.7          -0.7470588 4.410535e-04        4.889610e-03
#> hsa.miR.873.5p.11         -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.874.3p.5          -0.7500000 4.098820e-04        4.880632e-03
#> hsa.miR.874.3p.6          -0.7764706 2.021754e-04        3.484055e-03
#> hsa.miR.874.3p.7          -0.5558824 1.268109e-02        4.056686e-02
#> hsa.miR.874.3p.11         -0.5558824 1.268109e-02        4.056686e-02
#> hsa.miR.874.3p.15         -0.7470588 4.410535e-04        4.889610e-03
#> hsa.miR.874.3p.16         -0.6088235 6.157602e-03        2.342803e-02
#> hsa.miR.874.3p.17         -0.7705882 2.384112e-04        3.484055e-03
#> hsa.miR.9.5p.2            -0.7029412 1.193476e-03        7.993805e-03
#> hsa.miR.9.5p.3            -0.5647059 1.133160e-02        3.755784e-02
#> hsa.miR.9.5p.9            -0.6205882 5.156758e-03        2.046787e-02
#> hsa.miR.9.5p.12           -0.7470588 4.410535e-04        4.889610e-03
#> hsa.miR.96.5p.1           -0.6058824 6.430149e-03        2.432109e-02
#> hsa.miR.96.5p.2           -0.8029412 9.025840e-05        2.418173e-03
#> hsa.miR.96.5p.3           -0.6264706 4.706989e-03        1.991180e-02
#> hsa.miR.96.5p.5           -0.7588235 3.270151e-04        4.043668e-03
#> hsa.miR.96.5p.8           -0.7029412 1.193476e-03        7.993805e-03
#> hsa.miR.96.5p.9           -0.6970588 1.345295e-03        8.317543e-03
#> hsa.miR.96.5p.15          -0.6352941 4.091032e-03        1.814161e-02
#> hsa.miR.96.5p.18          -0.5882353 8.269861e-03        2.921715e-02
#> hsa.miR.96.5p.19          -0.6000000 7.003561e-03        2.618192e-02
#> hsa.miR.96.5p.20          -0.5823529 8.966367e-03        3.133355e-02
```

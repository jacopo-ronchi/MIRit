# basic differential expression works for edgeR

    Code
      mirnaDE(de_edgeR)
    Output
                               ID     logFC   AveExpr      P.Value    adj.P.Val
      hsa-let-7e-5p hsa-let-7e-5p  1.076945 13.575832 6.209083e-06 0.0002048997
      hsa-miR-1179   hsa-miR-1179 -1.953635  7.593888 6.810657e-04 0.0112375836
      hsa-miR-1         hsa-miR-1  2.335482  6.622710 7.439913e-03 0.0481280074

# basic differential expression works for DESeq2

    Code
      mirnaDE(de_DESeq2)
    Output
                             ID     logFC  AveExpr      P.Value   adj.P.Val
      hsa-miR-1179 hsa-miR-1179 -2.075882 76.06323 0.0002964329 0.006076874

# basic differential expression works for limma-voom

    Code
      mirnaDE(de_voom)
    Output
                               ID     logFC   AveExpr      P.Value    adj.P.Val
      hsa-let-7e-5p hsa-let-7e-5p  1.130509 13.437294 2.202550e-07 7.268414e-06
      hsa-miR-1179   hsa-miR-1179 -2.186194  6.745134 3.597985e-04 4.107650e-03

# basic differential expression works for limma

    Code
      mirnaDE(de_limma)
    Output
                             ID     logFC  AveExpr      P.Value  adj.P.Val
      hsa-miR-1275 hsa-miR-1275 -1.101810 3.164477 0.0002656976 0.01328488
      hsa-miR-1179 hsa-miR-1179 -2.403837 4.534627 0.0006524745 0.01631186


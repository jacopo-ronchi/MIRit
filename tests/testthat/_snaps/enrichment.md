# basic ORA enrichment works

    Code
      enrichmentResults(enr$downregulated)
    Output
                                              pathway         pval       padj overlap
      1             thyroid hormone metabolic process 1.155722e-05 0.02092748       6
      2  phenol-containing compound metabolic process 2.373480e-05 0.02092748       7
      3                  regulation of hormone levels 2.590035e-05 0.02092748      13
      4                     hormone metabolic process 4.190975e-05 0.02539731       9
      5                           response to hypoxia 1.025394e-04 0.04074207      11
      6          skeletal muscle cell differentiation 1.176545e-04 0.04074207       5
      7                    thyroid hormone generation 1.176545e-04 0.04074207       5
      8               response to endogenous stimulus 1.460546e-04 0.04234129      28
      9           response to decreased oxygen levels 1.668038e-04 0.04234129      11
      10           skeletal muscle tissue development 1.746753e-04 0.04234129       7
      11                    response to oxygen levels 2.264697e-04 0.04990568      11
         size overlapGenes
      1    12 PAX8, IY....
      2    19 TPO, PAX....
      3    68 FOXE1, T....
      4    35 TG, PAX8....
      5    57 APAF1, J....
      6    11 EGR1, EG....
      7    11 FOXE1, T....
      8   270 GJA1, UF....
      9    60 RYR2, AP....
      10   25 BCL2, AT....
      11   62 RYR2, EG....

# basic GSEA enrichment works

    Code
      enrichmentResults(enr)
    Output
                                            pathway         pval       padj   log2err
      1                   Thyroid hormone synthesis 5.572528e-05 0.01192521 0.5573322
      2 Protein processing in endoplasmic reticulum 2.271988e-04 0.02431027 0.5188481
      3              Coronavirus disease - COVID-19 6.085682e-04 0.04341120 0.4772708
                ES       NES size  leadingEdge
      1 -0.6739420 -2.017081   27 SLC26A4,....
      2 -0.5431852 -1.883585   55 BCL2, HY....
      3 -0.4496718 -1.676284   85 JUN, RPL....

# basic CAMERA enrichment works

    Code
      enrichmentResults(enr)
    Output
      [1] pathway   direction overlap   size      pval      padj     
      <0 rows> (or 0-length row.names)


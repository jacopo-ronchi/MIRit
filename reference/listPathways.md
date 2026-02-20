# List all the available biological pathways in KEGG, Reactome and WikiPathways

This function can be used to retrieve a list of valid biological
pathways present in KEGG, Reactome and WikiPathways.

## Usage

``` r
listPathways(organism, database)
```

## Arguments

- organism:

  The name of the organism under consideration. The different databases
  have different supported organisms. To see the list of supported
  organisms for a given database, use the
  [`supportedOrganisms()`](https://jacopo-ronchi.github.io/MIRit/reference/supportedOrganisms.md)
  function

- database:

  The name of the database to use. It must be one of: `KEGG`,
  `Reactome`, and `WikiPathways`

## Value

A `character` vector containing the pathway names present in the
specified database.

## Note

This function uses the `graphite` package to retrieve biological
pathways from KEGG, Reactome and WikiPathways.

## References

Sales, G., Calura, E., Cavalieri, D. et al. graphite - a Bioconductor
package to convert pathway topology to gene network. BMC Bioinformatics
13, 20 (2012), <https://doi.org/10.1186/1471-2105-13-20>.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# list the mouse pathways present in WikiPathways
listPathways("Mus musculus", "WikiPathways")
#>   [1] "Statin pathway"                                                                  
#>   [2] "Cholesterol biosynthesis"                                                        
#>   [3] "Selenium metabolism / selenoproteins"                                            
#>   [4] "TGF-beta signaling pathway"                                                      
#>   [5] "Hedgehog signaling pathway"                                                      
#>   [6] "Glucuronidation"                                                                 
#>   [7] "EBV LMP1 signaling"                                                              
#>   [8] "Estrogen signaling"                                                              
#>   [9] "Transcriptional activation by Nfe2l2 in response to phytochemicals"              
#>  [10] "Methylation"                                                                     
#>  [11] "EPO receptor signaling"                                                          
#>  [12] "Amino acid conjugation of benzoic acid"                                          
#>  [13] "Type II interferon signaling (IFNG)"                                             
#>  [14] "Apoptosis"                                                                       
#>  [15] "Nod-like receptor (NLR) signaling pathway"                                       
#>  [16] "Retinol metabolism"                                                              
#>  [17] "ErbB signaling pathway"                                                          
#>  [18] "Aflatoxin B1 metabolism"                                                         
#>  [19] "Mitochondrial gene expression"                                                   
#>  [20] "Estrogen metabolism"                                                             
#>  [21] "Polyol pathway"                                                                  
#>  [22] "SIDS susceptibility pathways"                                                    
#>  [23] "Endochondral ossification"                                                       
#>  [24] "Selenium micronutrient network"                                                  
#>  [25] "Folic acid network"                                                              
#>  [26] "Oxidation by cytochrome P450"                                                    
#>  [27] "Oxidative damage response"                                                       
#>  [28] "Dopaminergic neurogenesis"                                                       
#>  [29] "Regulation of cardiac hypertrophy by miR-208"                                    
#>  [30] "MicroRNAs in cardiomyocyte hypertrophy"                                          
#>  [31] "Glycolysis and gluconeogenesis"                                                  
#>  [32] "Iron homeostasis"                                                                
#>  [33] "Cytoplasmic ribosomal proteins"                                                  
#>  [34] "Glutathione metabolism"                                                          
#>  [35] "Apoptosis modulation by HSP70"                                                   
#>  [36] "Acetylcholine synthesis"                                                         
#>  [37] "Mechanisms associated with pluripotency"                                         
#>  [38] "One-carbon metabolism and related pathways"                                      
#>  [39] "Kennedy pathway"                                                                 
#>  [40] "Heme biosynthesis"                                                               
#>  [41] "GPCRs, class A rhodopsin-like"                                                   
#>  [42] "Hepatocyte growth factor receptor signaling"                                     
#>  [43] "Splicing factor NOVA regulated synaptic proteins"                                
#>  [44] "Complement activation, classical pathway"                                        
#>  [45] "Ptf1a related regulatory pathway"                                                
#>  [46] "Hypertrophy model"                                                               
#>  [47] "Heart development"                                                               
#>  [48] "Neural crest differentiation"                                                    
#>  [49] "Alzheimer's disease"                                                             
#>  [50] "Serotonin receptor 2 and STAT3 signaling"                                        
#>  [51] "SREBF and miR33 in cholesterol and lipid homeostasis"                            
#>  [52] "Serotonin and anxiety-related events"                                            
#>  [53] "Serotonin and anxiety"                                                           
#>  [54] "BDNF pathway"                                                                    
#>  [55] "Striated muscle contraction"                                                     
#>  [56] "Purine metabolism"                                                               
#>  [57] "Chemokine signaling pathway"                                                     
#>  [58] "PPAR signaling pathway"                                                          
#>  [59] "Fatty acid oxidation"                                                            
#>  [60] "G protein signaling pathways"                                                    
#>  [61] "miRNAs and TFs in iPS Cell Generation"                                           
#>  [62] "Osteoblast signaling"                                                            
#>  [63] "Spinal cord injury"                                                              
#>  [64] "TNF-alpha NF-kB signaling pathway"                                               
#>  [65] "Mapk cascade"                                                                    
#>  [66] "Primary focal segmental glomerulosclerosis (FSGS)"                               
#>  [67] "Focal adhesion: PI3K-Akt-mTOR signaling pathway"                                 
#>  [68] "Gene regulatory network modelling somitogenesis"                                 
#>  [69] "White fat cell differentiation"                                                  
#>  [70] "Notch signaling pathway"                                                         
#>  [71] "Electron transport chain"                                                        
#>  [72] "G13 signaling pathway"                                                           
#>  [73] "Translation factors"                                                             
#>  [74] "Glycogen metabolism"                                                             
#>  [75] "Eicosanoid synthesis"                                                            
#>  [76] "Fatty acid omega-oxidation"                                                      
#>  [77] "ESC pluripotency pathways"                                                       
#>  [78] "p38 Mapk signaling pathway"                                                      
#>  [79] "ApoE and miR-146 in inflammation and atherosclerosis"                            
#>  [80] "Tyrobp causal network in microglia"                                              
#>  [81] "Microglia pathogen phagocytosis pathway"                                         
#>  [82] "Lung fibrosis"                                                                   
#>  [83] "Parkinson's disease"                                                             
#>  [84] "Ectodysplasin A signaling in hair follicle development"                          
#>  [85] "Novel Jun-Dmp1 pathway"                                                          
#>  [86] "BMP signaling pathway in eyelid development"                                     
#>  [87] "Hfe effect on hepcidin production"                                               
#>  [88] "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling"  
#>  [89] "IL-1 signaling pathway"                                                          
#>  [90] "Prostaglandin synthesis and regulation"                                          
#>  [91] "Myometrial relaxation and contraction pathways"                                  
#>  [92] "Wnt signaling in kidney disease"                                                 
#>  [93] "Robo4 and VEGF signaling pathways crosstalk"                                     
#>  [94] "ACE inhibitor pathway"                                                           
#>  [95] "miR-127 in mesendoderm differentiation"                                          
#>  [96] "Wnt signaling"                                                                   
#>  [97] "Oxidative stress response"                                                       
#>  [98] "G1 to S cell cycle control"                                                      
#>  [99] "Distal convoluted tubule 1 (DCT1) cell"                                          
#> [100] "Ethanol metabolism resulting in production of ROS by CYP2E1"                     
#> [101] "Nuclear receptors in lipid metabolism and toxicity"                              
#> [102] "Eicosanoid lipid synthesis map"                                                  
#> [103] "TCA cycle"                                                                       
#> [104] "Sphingolipid metabolism overview"                                                
#> [105] "Glycerolipids and glycerophospholipids"                                          
#> [106] "Cholesterol metabolism with Bloch and Kandutsch-Russell pathways"                
#> [107] "Eicosanoid metabolism via cyclooxygenases (COX)"                                 
#> [108] "Eicosanoid metabolism via lipoxygenases (LOX)"                                   
#> [109] "Eicosanoid metabolism via cytochrome P450 monooxygenases"                        
#> [110] "One-carbon metabolism"                                                           
#> [111] "Omega-3 / omega-6 fatty acid synthesis"                                          
#> [112] "Omega-9 fatty acid synthesis"                                                    
#> [113] "Oxidative stress and redox pathway"                                              
#> [114] "Circulating monocytes and cardiac macrophages in diastolic dysfunction"          
#> [115] "Complement and coagulation cascades"                                             
#> [116] "Elongation of (very) long chain fatty acids"                                     
#> [117] "Osteoclast signaling"                                                            
#> [118] "Inflammatory response pathway"                                                   
#> [119] "Blood clotting cascade"                                                          
#> [120] "Lipids measured in liver metastasis from breast cancer"                          
#> [121] "Sphingolipid metabolism (integrated pathway)"                                    
#> [122] "Regulation of Pgc1a expression by a Gsk3b-Tfeb signaling axis in skeletal muscle"
#> [123] "GDNF/RET signaling axis"                                                         
#> [124] "Peroxiredoxin 2 induced ovarian failure"                                         
#> [125] "Mapk signaling pathway"                                                          
#> [126] "Deregulation of renin-angiotensin system by SARS-CoV infection"                  
#> [127] "Hypoxia-dependent self-renewal of myoblasts"                                     
#> [128] "Hypoxia-dependent proliferation of myoblasts"                                    
#> [129] "Hypoxia-dependent differentiation of myoblasts"                                  
#> [130] "Na/K-ATPase/Src signaling"                                                       
#> [131] "Burn wound healing"                                                              
#> [132] "Fibrin complement receptor 3 signaling pathway"                                  
#> [133] "Oxylipins pathways"                                                              
#> [134] "Proteasome degradation"                                                          
#> [135] "Biogenic amine synthesis"                                                        
#> [136] "Regulation of actin cytoskeleton"                                                
#> [137] "Lac-Phe pathway"                                                                 
#> [138] "Comprehensive IL-17A signaling"                                                  
#> [139] "Dravet syndrome: Scn1a-A1783V point mutation model"                              
#> [140] "Globo series sphingolipid metabolism"                                            
#> [141] "Synthesis and degradation of ketone bodies"                                      
#> [142] "Exercise-induced circadian regulation"                                           
#> [143] "Steroid biosynthesis"                                                            
#> [144] "Calcium regulation in cardiac cells"                                             
#> [145] "Signal transduction of S1P receptor"                                             
#> [146] "FAS pathway and stress induction of HSP regulation"                              
#> [147] "Leptin-insulin signaling overlap"                                                
#> [148] "Integrin-mediated cell adhesion"                                                 
#> [149] "miR-1 in cardiac development"                                                    
#> [150] "Pentose phosphate pathway"                                                       
#> [151] "Insulin signaling"                                                               
#> [152] "Amino acid metabolism"                                                           
#> [153] "Leptin and adiponectin"                                                          
#> [154] "Wnt signaling pathway and pluripotency"                                          
#> [155] "Glutathione and one-carbon metabolism"                                           
#> [156] "Focal adhesion"                                                                  
#> [157] "Toll-like receptor signaling"                                                    
#> [158] "Oxidative phosphorylation"                                                       
#> [159] "Arachidonate epoxygenase / epoxide hydrolase"                                    
#> [160] "Metapathway biotransformation"                                                   
#> [161] "Fatty acid beta-oxidation"                                                       
#> [162] "Fatty acid biosynthesis"                                                         
#> [163] "Tryptophan metabolism"                                                           
#> [164] "Nucleotide GPCRs"                                                                
#> [165] "GPCRs, small ligand"                                                             
#> [166] "Monoamine GPCRs"                                                                 
```

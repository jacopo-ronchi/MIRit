# Get the list of supported organisms for a given database

This function provides the list of supported organisms for different
databases, namely Gene Ontology (GO), Kyoto Encyclopedia of Genes and
Genomes (KEGG), MsigDB, WikiPathways, Reactome, Enrichr, Disease
Ontology (DO), Network of Cancer Genes (NCG), DisGeNET, and COVID19.

## Usage

``` r
supportedOrganisms(database)
```

## Arguments

- database:

  The database name. It must be one of: `GO`, `KEGG`, `MsigDB`,
  `WikiPathways`, `Reactome`, `Enrichr`, `DO`, `NCG`, `DisGeNET`,
  `COVID19`

## Value

A `character` vector listing all the supported organisms for the
database specified by the user.

## Note

To perform the functional enrichment of genes, MIRit uses the `geneset`
R package to download gene sets from the above mentioned databases.

## References

Liu, Y., Li, G. Empowering biologists to decode omics data: the Genekitr
R package and web server. BMC Bioinformatics 24, 214 (2023).
<https://doi.org/10.1186/s12859-023-05342-9>.

## Author

Jacopo Ronchi, <jacopo.ronchi@unimib.it>

## Examples

``` r
# get the supported organisms for GO database
supportedOrganisms("GO")
#>   [1] "Amborella trichopoda"               "Anolis carolinensis"               
#>   [3] "Anopheles gambiae"                  "Aquifex aeolicus"                  
#>   [5] "Arabidopsis thaliana"               "Ashbya gossypii"                   
#>   [7] "Bacillus cereus"                    "Bacillus subtilis"                 
#>   [9] "Bacteroides thetaiotaomicron"       "Batrachochytrium dendrobatidis"    
#>  [11] "Bos taurus"                         "Brachypodium distachyon"           
#>  [13] "Bradyrhizobium diazoefficiens"      "Branchiostoma floridae"            
#>  [15] "Brassica napus"                     "Brassica rapa subsp. pekinensis"   
#>  [17] "Caenorhabditis briggsae"            "Caenorhabditis elegans"            
#>  [19] "Candida albicans"                   "Canis lupus familiaris"            
#>  [21] "Capsicum annuum"                    "Chlamydia trachomatis"             
#>  [23] "Chlamydomonas reinhardtii"          "Chloroflexus aurantiacus"          
#>  [25] "Ciona intestinalis"                 "Citrus sinensis"                   
#>  [27] "Clostridium botulinum"              "Coxiella burnetii"                 
#>  [29] "Cryptococcus neoformans"            "Cucumis sativus"                   
#>  [31] "Danio rerio"                        "Daphnia pulex"                     
#>  [33] "Deinococcus radiodurans"            "Dictyoglomus turgidum"             
#>  [35] "Dictyostelium discoideum"           "Dictyostelium purpureum"           
#>  [37] "Drosophila melanogaster"            "Emericella nidulans"               
#>  [39] "Entamoeba histolytica"              "Equus caballus"                    
#>  [41] "Erythranthe guttata"                "Escherichia coli"                  
#>  [43] "Eucalyptus grandis"                 "Felis catus"                       
#>  [45] "Fusobacterium nucleatum"            "Gallus gallus"                     
#>  [47] "Geobacter sulfurreducens"           "Giardia intestinalis"              
#>  [49] "Gloeobacter violaceus"              "Glycine max"                       
#>  [51] "Gorilla gorilla gorilla"            "Gossypium hirsutum"                
#>  [53] "Haemophilus influenzae"             "Halobacterium salinarum"           
#>  [55] "Helianthus annuus"                  "Helicobacter pylori"               
#>  [57] "helobdella robusta"                 "Homo sapiens"                      
#>  [59] "Hordeum vulgare subsp. vulgare"     "Ixodes scapularis"                 
#>  [61] "Juglans regia"                      "Klebsormidium nitens"              
#>  [63] "Korarchaeum cryptofilum"            "Lactuca sativa"                    
#>  [65] "Leishmania major"                   "lepisosteus oculatus"              
#>  [67] "Leptospira interrogans"             "Listeria monocytogenes"            
#>  [69] "Macaca mulatta"                     "Manihot esculenta"                 
#>  [71] "Marchantia polymorpha"              "Medicago truncatula"               
#>  [73] "Methanocaldococcus jannaschii"      "Methanosarcina acetivorans"        
#>  [75] "Monodelphis domestica"              "Monosiga brevicollis"              
#>  [77] "Mus musculus"                       "Musa acuminata subsp. malaccensis" 
#>  [79] "Mycobacterium tuberculosis"         "mycoplasma genitalium"             
#>  [81] "Neisseria meningitidis serogroup b" "Nelumbo nucifera"                  
#>  [83] "Nematostella vectensis"             "Neosartorya fumigata"              
#>  [85] "Neurospora crassa"                  "Nicotiana tabacum"                 
#>  [87] "Nitrosopumilus maritimus"           "Ornithorhynchus anatinus"          
#>  [89] "Oryza sativa"                       "Oryzias latipes"                   
#>  [91] "Pan troglodytes"                    "Paramecium tetraurelia"            
#>  [93] "Phaeosphaeria nodorum"              "Physcomitrella patens"             
#>  [95] "Phytophthora ramorum"               "Plasmodium falciparum"             
#>  [97] "Populus trichocarpa"                "Pristionchus pacificus"            
#>  [99] "Prunus persica"                     "Pseudomonas aeruginosa"            
#> [101] "Puccinia graminis"                  "Pyrobaculum aerophilum"            
#> [103] "Rattus norvegicus"                  "Rhodopirellula baltica"            
#> [105] "Ricinus communis"                   "Saccharomyces cerevisiae"          
#> [107] "Salmonella typhimurium"             "Schizosaccharomyces japonicus"     
#> [109] "Schizosaccharomyces pombe"          "Sclerotinia sclerotiorum"          
#> [111] "Selaginella moellendorffii"         "Setaria italica"                   
#> [113] "Shewanella oneidensis"              "Solanum lycopersicum"              
#> [115] "Solanum tuberosum"                  "Sorghum bicolor"                   
#> [117] "Spinacia oleracea"                  "Staphylococcus aureus"             
#> [119] "Streptococcus pneumoniae"           "Streptomyces coelicolor"           
#> [121] "Strongylocentrotus purpuratus"      "Sulfolobus solfataricus"           
#> [123] "Sus scrofa"                         "Synechocystis"                     
#> [125] "Thalassiosira pseudonana"           "Theobroma cacao"                   
#> [127] "Thermococcus kodakaraensis"         "Thermodesulfovibrio yellowstonii"  
#> [129] "Thermotoga maritima"                "Tribolium castaneum"               
#> [131] "Trichomonas vaginalis"              "Trichoplax adhaerens"              
#> [133] "Triticum aestivum"                  "Trypanosoma brucei"                
#> [135] "Ustilago maydis"                    "Vibrio cholerae"                   
#> [137] "Vitis vinifera"                     "Xanthomonas campestris"            
#> [139] "Xenopus tropicalis"                 "Yarrowia lipolytica"               
#> [141] "Yersinia pestis"                    "Zea mays"                          
#> [143] "Zostera marina"                    

# get the supported organisms for Reactome
supportedOrganisms("Reactome")
#>  [1] "Bos taurus"               "Caenorhabditis elegans"  
#>  [3] "Danio rerio"              "Drosophila melanogaster" 
#>  [5] "Gallus gallus"            "Homo sapiens"            
#>  [7] "Mus musculus"             "Rattus norvegicus"       
#>  [9] "Saccharomyces cerevisiae" "Sus scrofa"              
#> [11] "Xenopus tropicalis"      
```

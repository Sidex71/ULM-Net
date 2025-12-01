# adding cell assignments to seurat object

'AddMtaObject()' a function to add the predicted cell labels to the
metadata of scRNAseq

## Usage

``` r
AddMetaObject(seurat_obj, cell_class_df)
```

## Arguments

- seurat_obj:

  a prepossessed Seurat object storing the scRNAseq data

- cell_class_df:

  a data frame of cell barcodes and cell type assignments, ideally
  obtained as the output from the GetCellAssignments() function.

## Value

a new Seurat object with the updated metadata containing predicted cell
labels in the "celltype_ulm" column

## Examples

``` r
data(int_multData)
data(int_signature)
my_scores <- GetCellScores(seurat_obj = int_multData[,1:1000], signatures = int_signature, assay = 'RNA', layer = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
new_obj <- AddMetaObject(seurat_obj = int_multData[,1:1000], cell_class_df = my_ass)
head(new_obj$celltype_ulm, 20)
#>                      AAACCCACAAATCGTC                      AAACCCACACAACGCC 
#>                    "Progenitor.early"                   "Progenitor.late.2" 
#>                      AAACCCAGTATCTCGA                      AAACCCAGTCATGCAT 
#>                              "Goblet"                              "Goblet" 
#>                      AAACGAAAGGTAGTCG                      AAACGAACACCAACAT 
#>                              "Goblet"                       "Goblet_Paneth" 
#>                      AAACGAAGTATCGAAA                      AAACGAATCCGTAGGC 
#>                                    NA        "Enterocyte_Progenitor.late.1" 
#>                      AAACGCTCAACGGGTA                      AAACGCTCAGCCTACG 
#>        "Enterocyte_Progenitor.late.1" "Progenitor.early_Transit.amplifying" 
#>                      AAACGCTGTATGAGAT                      AAACGCTTCATGGAGG 
#>                    "Progenitor.early"                              "Goblet" 
#>                      AAACGCTTCTTCTAAC                      AAAGAACAGACATCAA 
#>                                "Tuft"                              "Goblet" 
#>                      AAAGAACAGGCTTCCG                      AAAGAACCAGAGACTG 
#> "Progenitor.late.1_Progenitor.late.2"                              "Goblet" 
#>                      AAAGAACCAGCAGAAC                      AAAGAACCAGCGACCT 
#>                                    NA                              "Goblet" 
#>                      AAAGAACGTCCACAGC                      AAAGGATGTCATCAGT 
#>                   "Progenitor.late.1"                              "Goblet" 
```

# obtaining multiplet data

'GetMultiplet()' a function to extract predicted multiplets

## Usage

``` r
GetMultiplet(seurat_obj, minCells = 2)
```

## Arguments

- seurat_obj:

  a Seurat object with the metadata containing predicted cell labels in
  the "celltype_ulm" column and the number of cells in the "count_ulm"
  column. Ideally the output from the AddMetObject() function

- minCells:

  a numeric value specifying the minimum number of cells. Default is 2
  to include doublets and/or higher order multiplets

## Value

a list containing a Seurat object of multiplets and multiplet
distribution summary

## Examples

``` r
data(int_multData)
data(int_signature)
my_scores <- GetCellScores(seurat_obj = int_multData[,1:1000], signatures = int_signature, assay = 'RNA', layer = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
new_obj <- AddMetaObject(seurat_obj = int_multData[,1:1000], cell_class_df = my_ass)
my_mult <- GetMultiplet(seurat_obj = new_obj)
#> Warning: Removing 223 cells missing data for vars requested
my_mult
#> $multSummary
#>                                            multipletType frequency
#> 1  Enterocyte_Goblet_Progenitor.late.1_Progenitor.late.2         1
#> 2                    Enterocyte_Paneth_Progenitor.late.1         1
#> 3                           Enterocyte_Progenitor.late.1        69
#> 4         Enterocyte_Progenitor.late.1_Progenitor.late.2        15
#> 5                      Enterocyte_Progenitor.late.1_Tuft         1
#> 6                     Enteroendocrine_Transit.amplifying         1
#> 7                                          Goblet_Paneth        56
#> 8                                Goblet_Progenitor.early        24
#> 9             Goblet_Progenitor.early_Transit.amplifying         2
#> 10                              Goblet_Progenitor.late.1         1
#> 11            Goblet_Progenitor.late.1_Progenitor.late.2         1
#> 12                              Goblet_Progenitor.late.2         2
#> 13                             Goblet_Transit.amplifying         9
#> 14                               Paneth_Progenitor.early         5
#> 15                              Paneth_Progenitor.late.1         1
#> 16            Paneth_Progenitor.late.1_Progenitor.late.2         1
#> 17                                           Paneth_Stem         2
#> 18                             Paneth_Transit.amplifying         4
#> 19                    Progenitor.early_Progenitor.late.1         1
#> 20                   Progenitor.early_Transit.amplifying        21
#> 21                   Progenitor.late.1_Progenitor.late.2        57
#> 22              Progenitor.late.1_Progenitor.late.2_Tuft         2
#> 
#> $multObj
#> An object of class Seurat 
#> 15615 features across 277 samples within 1 assay 
#> Active assay: RNA (15615 features, 2000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
#> 


```

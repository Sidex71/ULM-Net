# filtering multiplets

'FilterMultiplet()' a function to filter multiplets of a predefined
frequency

## Usage

``` r
FilterMultiplet(seurat_obj, minCells = 2, minFreq = 10)
```

## Arguments

- seurat_obj:

  a Seurat object with the metadata containing predicted cell labels in
  the "celltype_ulm" column and the number of cells in the "count_ulm"
  column. Ideally the output from the AddMetObject() function or the
  multiplet Seurat object from the GetMultiplet() function

- minCells:

  a numeric value specifying the minimum number of cells. Default is 2
  to include doublets and/or higher order multiplets

- minFreq:

  a numeric value specifying the minimum frequency of a multiplet type
  for it to be retained. Default is 10.

## Value

a filtered list containing a Seurat object of multiplets and a dataframe
of multiplet distribution summary.

## Examples

``` r
data(int_multData)
data(int_signature)
my_scores <- GetCellScores(seurat_obj = int_multData[,1:1000], signatures = int_signature, assay = 'RNA', layer = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
new_obj <- AddMetaObject(seurat_obj = int_multData[,1:1000], cell_class_df = my_ass)
my_mult <- GetMultiplet(seurat_obj = new_obj)
#> Warning: Removing 223 cells missing data for vars requested
my_mult_filt <- FilterMultiplet(seurat_obj = new_obj)
#> Warning: Removing 223 cells missing data for vars requested
my_mult_filt
#> $multSummaryFilt
#>                                     multipletType frequency
#> 3                    Enterocyte_Progenitor.late.1        69
#> 4  Enterocyte_Progenitor.late.1_Progenitor.late.2        15
#> 7                                   Goblet_Paneth        56
#> 8                         Goblet_Progenitor.early        24
#> 20            Progenitor.early_Transit.amplifying        21
#> 21            Progenitor.late.1_Progenitor.late.2        57
#> 
#> $multObjFilt
#> An object of class Seurat 
#> 15615 features across 242 samples within 1 assay 
#> Active assay: RNA (15615 features, 2000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
#> 
```

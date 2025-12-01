# getting final cell assignments

'GetCellAssignments()' a function that assigns cell type labels to each
barcode based on the signature scores

## Usage

``` r
GetCellAssignments(score_data, p_val = 0.05, cut_off = 1)
```

## Arguments

- score_data:

  a data frame of cell barcodes and gene signature scores, ideally
  obtained as the output from the GetCellScores() function.

- p_val:

  a numeric value specifying the p value cut-off to filter significant
  signature scores (default: 0.05)

- cut_off:

  a numeric value specifying the cut-off for signature scores (default:
  1)

## Value

a data frame of cell barcodes and cell type assignments. A barcode may
be assigned a single or multi cell type assignment depending on
signature enrichment scores.

## Examples

``` r
data(int_singData)
data(int_signature)
my_scores <- GetCellScores(seurat_obj = int_singData[,1:1000], signatures = int_signature, assay = 'RNA', layer = 'data')
my_ass <- GetCellAssignments(score_data = my_scores)
head(my_ass)
#>            barcode statistic count_ulm                        celltype_ulm
#> 1 AAAGGTATCGGCTTCT       ulm         2        Enterocyte_Progenitor.late.1
#> 2 AACAGGGAGAAGTCAT       ulm         2        Enterocyte_Progenitor.late.1
#> 3 AACCAACCATCAGCAT       ulm         2        Enterocyte_Progenitor.late.1
#> 4 AACCTGAGTACCGTCG       ulm         3 Enterocyte_Paneth_Progenitor.late.1
#> 5 AACGAAACAAGCGCAA       ulm         2        Enterocyte_Progenitor.late.1
#> 6 AACGAAACAGGTCCCA       ulm         2        Enterocyte_Progenitor.late.1
#>     avg_pvalue avg_score
#> 1 2.523098e-13  8.319202
#> 2 2.424298e-15  8.761571
#> 3 9.232242e-16  8.127915
#> 4 8.578408e-06  7.990285
#> 5 2.859176e-09  5.994553
#> 6 1.859488e-17  8.907267
```

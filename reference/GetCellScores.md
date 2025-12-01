# scoring cells for gene signatures

'GetCellScores()' scores each cell in the scRNAseq data for cell
type-specific gene signatures

## Usage

``` r
GetCellScores(
  seurat_obj,
  signatures,
  assay = "RNA",
  slot = NULL,
  layer = "data"
)
```

## Arguments

- seurat_obj:

  a preprocessed Seurat object storing the scRNAseq data

- signatures:

  a data frame of cell type signatures with cell types in the source
  column, genes in the target column, and weights in the mor column.
  Ideally the output from the GetSignature() function.

- assay:

  a character specifying assay in the Seurat object containing count
  matrix.

- layer:

  a character specifying the layer to draw counts from (Seurat v3/v4).
  This can be "counts" for raw counts, "data" for normalized counts, or
  "scaled" for scaled counts. Default is NULL.

## Value

a data frame of cell barcodes and gene signature scores

## Examples

``` r
data(int_singData)
data(int_signature)
my_scores <- GetCellScores(seurat_obj = int_singData[,1:1000],
                           signatures = int_signature,
                           assay = 'RNA',
                           layer = 'data')
head(my_scores)
#> # A tibble: 6 × 5
#>   barcode          celltype   score   p_value statistic
#>   <chr>            <chr>      <dbl>     <dbl> <chr>    
#> 1 AAACGAAAGAGGTCGT Enterocyte -3.20 0.00136   ulm      
#> 2 AAACGAAAGCCTCCAG Enterocyte -3.30 0.000974  ulm      
#> 3 AAACGAAAGTTTGGCT Enterocyte -3.88 0.000106  ulm      
#> 4 AAACGAACACAAGCTT Enterocyte -4.31 0.0000164 ulm      
#> 5 AAACGAAGTCGCTTGG Enterocyte -1.22 0.221     ulm      
#> 6 AAACGAAGTCTTACTT Enterocyte -4.01 0.0000622 ulm      
```

# generating gene signature

'GetSignature()' generates cell type-specific gene signatures from
scRNAseq data

## Usage

``` r
GetSignature(seurat_obj, ident_col = NULL, n = 100, p_val = 0.05)
```

## Arguments

- seurat_obj:

  a prepossessed Seurat object storing the scRNAseq data

- ident_col:

  a column in the Seurat object metadata which is a character vector
  storing cell names or labels. If not specified, the default ident of
  the Seurat object will be used.

- n:

  a numeric value specifying the number of genes to be used for each
  cell signature. Default is 100 genes per cell type.

- p_val:

  a numeric value specifying the adjusted p-value cut-off

## Value

a dataframe of cell type signatures.

## Examples

``` r
data(int_singData)
int_sig <- GetSignature(seurat_obj = int_singData[,1:1000], ident_col = int_singData$Cell_Type)
#> using the specified seurat ident to generate signatures
#> Calculating cluster Progenitor early
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#> Calculating cluster Progenitor late-1
#> Calculating cluster Transit amplifying
#> Calculating cluster Progenitor late-2
#> Calculating cluster Goblet
#> Calculating cluster Stem
#> Calculating cluster Enterocyte
#> Calculating cluster Paneth
#> Calculating cluster Enteroendocrine
#> Calculating cluster Tuft
head(int_sig)
#> # A tibble: 6 × 3
#> # Groups:   source [1]
#>   source           target          mor
#>   <chr>            <chr>         <dbl>
#> 1 Progenitor early C330021F23Rik     1
#> 2 Progenitor early Cdc25c            1
#> 3 Progenitor early Knstrn            1
#> 4 Progenitor early Ccnb2             1
#> 5 Progenitor early Cdkn3             1
#> 6 Progenitor early Cenpa             1
```

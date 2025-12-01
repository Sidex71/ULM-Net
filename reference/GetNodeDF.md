# getting nodes and edges

'GetNodeDF()' a function to generate nodes and edges from multiplet data

## Usage

``` r
GetNodeDF(mat)
```

## Arguments

- mat:

  a data frame of multiplets types and their frequencies, ideally the
  summary data frame obtained from the FilterMultiplet() or
  GetMultiplet() function.

## Value

a data frame of cell-cell edges and nodes and the corresponding
frequencies

## Examples

``` r
mat <- data.frame('multilpetType' = c(paste('A', 'B', 'C', sep = '_'),
                                      paste('A', 'B', sep = '_'),
                                      paste('B', 'C', 'A', 'E', sep = '_'),
                                      paste('D', 'C', 'E', sep = '_')),
                  frequency = rep(50, 4))

mat
#>   multilpetType frequency
#> 1         A_B_C        50
#> 2           A_B        50
#> 3       B_C_A_E        50
#> 4         D_C_E        50
GetNodeDF(mat)
#>    Cell1 Cell2 n_cells
#> 1      B     A      50
#> 2      C     A      50
#> 3      A     B     100
#> 4      A     C      50
#> 5      B     C     100
#> 6      D     C      50
#> 7      A     E      50
#> 8      B     E      50
#> 9      C     E     100
#> 10     D     E      50
```

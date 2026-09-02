# List available target-strength models

List available target-strength models

## Usage

``` r
available_models()
```

## Value

A data frame describing currently available built-in and user-registered
target-strength models.

## Examples

``` r
head(available_models())
#>                   model        slot  source persistent aliases
#> bbfm               bbfm        BBFM builtin      FALSE        
#> bcms               bcms        BCMS builtin      FALSE        
#> calibration calibration calibration builtin      FALSE   soems
#> dwba               dwba        DWBA builtin      FALSE        
#> dwba_curved dwba_curved DWBA_curved builtin      FALSE        
#> ecms               ecms        ECMS builtin      FALSE        
```

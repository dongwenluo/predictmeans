# Multiple Comparisons Based on the Confidence Intervals

This function produces letter representations for a multiple comparison
test by analyzing the confidence intervals associated with the mean
values of different treatments. In particular, if the confidence
intervals of two treatments overlap, it indicates that there is no
significant difference between them. Conversely, if the confidence
intervals do not overlap, it indicates that the treatments are
significantly different from each other.

## Usage

``` r
ci_mcp(LL, UL, trt_n=NULL)
```

## Arguments

- LL:

  Lower limits of treatments' confidence interval.

- UL:

  Upper limits of treatments' confidence interval.

- trt_n:

  Treatments' names.

## References

Vanessa, C. (05 October 2022), *Confidence tricks: the 83.4% confidence
interval for comparing means*,
https://vsni.co.uk/blogs/confidence_trick.

## Author

Dongwen Luo, Siva Ganesh and John Koolaard

## Examples

``` r
  library(predictmeans)
  ci_mcp(LL=c(68.2566,  87.7566, 103.0899, 112.2566), UL=c(90.5212, 110.0212, 125.3545, 134.5212))
#>     1     2     3     4 
#> "A  " "AB " " BC" "  C" 
  
  data("Oats", package="nlme")
  Oats$nitro <- factor(Oats$nitro)
  fm <- lme(yield ~ nitro*Variety, random=~1|Block/Variety, data=Oats)
# fm <- lmer(yield ~ nitro*Variety+(1|Block/Variety), data=Oats)
  predictmeans(fm, "nitro", adj="BH", plot=FALSE)$mean_table
#> $Table
#>   nitro     Mean      SE   Df  LL(95%)  UL(95%) LetterGrp
#> 1     0  79.3889 7.17468 6.79  62.3165  96.4613      A   
#> 2   0.2  98.8889 7.17468 6.79  81.8165 115.9613       B  
#> 3   0.4 114.2222 7.17468 6.79  97.1498 131.2946        C 
#> 4   0.6 123.3889 7.17468 6.79 106.3165 140.4613         D
#> 
#> $Note
#> [1] "Letter-based representation of pairwise comparisons at significant level '0.05'"
#> 
  predictmeans(fm, "nitro", pair=TRUE, level=0.166, letterCI = TRUE, plot=FALSE)$mean_table
#> $Table
#>   nitro     Mean      SE   Df LL(83.4%) UL(83.4%) LetterGrp
#> 1     0  79.3889 7.17468 6.79   68.2563   90.5215       A  
#> 2   0.2  98.8889 7.17468 6.79   87.7563  110.0215       AB 
#> 3   0.4 114.2222 7.17468 6.79  103.0897  125.3548        BC
#> 4   0.6 123.3889 7.17468 6.79  112.2563  134.5215         C
#> 
#> $Note
#> [1] "Letter-based representation of pairwise comparisons by overlap of the 83.4% CIs (about significant level '0.05' for two uncorrelated means with equal SEs)"
#> 
```

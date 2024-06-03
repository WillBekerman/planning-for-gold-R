# planning.for.gold

R package to run the Sens-Val procedure, introduced by Bekerman et al. (2024), on pair-matched data from an observational study. The method identifies and tests hypotheses from a split sample that are robust to potential unmeasured biases.

# How to Install

You can install this through use devtools:

```r
devtools::install_github("WillBekerman/planning-for-gold-R", build_vignettes = TRUE)
```

# Example Usage

To see a real-data example of how to use the main function:

```r
vignette("Example", package = "planning.for.gold")
```

# Overview

For the full list of exported data and functions:

```{r}
library("planning.for.gold")
ls("package:planning.for.gold")
```

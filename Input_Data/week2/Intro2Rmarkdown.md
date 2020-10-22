Intro2Rmarkdown
================
Samantha Chen
10/21/2020

# title

## title

### title

basic text,

*italicize* single asterisk

**bold** double asterisks

1.  make
2.  numbered
3.  lists

we can make: - unnumbered - lists

[link here](https://google.com)

# Code Chunk

to add code, make a chunk by clicking (Insert and then R)

``` r
library(tidyverse)
```

    ## -- Attaching packages ----------

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts -------------------
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
# install.packages("RColorBrewer")
library(RColorBrewer)
```

    ## Warning: package 'RColorBrewer' was built under R version 4.0.3

or by using shortcuts, on mac (cmd +alt + i) and on pc (ctrl + alt + i)

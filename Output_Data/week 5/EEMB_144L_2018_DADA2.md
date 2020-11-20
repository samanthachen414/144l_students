EEMB\_144L\_2018\_DADA2
================
Samantha Chen
11/18/2020

``` r
library(tidyverse)
library(dada2)
library(ShortRead)
```

``` r
path <- "~/GitHub EEMB 144L/144l_students/Input_Data/week5/EEMB144L_2018_fastq/"

fnFs <- list.files(path, pattern = "_R1_001.fastq", full.names = TRUE)
fnRs <- list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)
```

``` r
FWD = "GTGYCAGCMGCCGCGGTAA"
REV = "GGACTACNVGGGTWTCTAAT"

all0rients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

FWD.orients <- all0rients(FWD)
REV.orients <- all0rients(REV)
```

``` r
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
```

    ##                  Forward Complement Reverse RevComp
    ## FWD.ForwardReads       0          0       0       0
    ## FWD.ReverseReads       0          0       0     283
    ## REV.ForwardReads       0          0       0    1195
    ## REV.ReverseReads       0          0       0       0

## Forward Reads

``` r
plotQualityProfile(fnFs[1:12])
```

![](EEMB_144L_2018_DADA2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Reverse Reads

``` r
plotQualityProfile(fnRs[1:12])
```

![](EEMB_144L_2018_DADA2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
\# Filtering and Trimming

``` r
sample.names <- sapply(strsplit(basename(fnFs), "_L"), '[', 1)
sample.names
```

    ##  [1] "144_A0_S6"  "144_A4_S7"  "144_A8_S8"  "144_B0_S9"  "144_B4_S10"
    ##  [6] "144_B8_S11" "144_C0_S12" "144_C4_S13" "144_C8_S14" "144_D0_S15"
    ## [11] "144_D4_S16" "144_D8_S17" "144_E0_S18" "144_E4_S19" "144_E8_S20"
    ## [16] "144_F0_S21" "144_F4_S22" "144_F8_S23" "144_G0_S24" "144_G4_S25"
    ## [21] "144_G8_S26" "144_H0_S27" "144_H4_S28" "144_H8_S29"

``` r
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(250,150), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE)

out
```

    ##                              reads.in reads.out
    ## 144_A0_S6_L001_R1_001.fastq     79521     66767
    ## 144_A4_S7_L001_R1_001.fastq     36911     31704
    ## 144_A8_S8_L001_R1_001.fastq     50351     42320
    ## 144_B0_S9_L001_R1_001.fastq     14331      3011
    ## 144_B4_S10_L001_R1_001.fastq    50889     42935
    ## 144_B8_S11_L001_R1_001.fastq    62383     52855
    ## 144_C0_S12_L001_R1_001.fastq    33962     23605
    ## 144_C4_S13_L001_R1_001.fastq    69734     59352
    ## 144_C8_S14_L001_R1_001.fastq    60793     52093
    ## 144_D0_S15_L001_R1_001.fastq    62933     55336
    ## 144_D4_S16_L001_R1_001.fastq    49383     42224
    ## 144_D8_S17_L001_R1_001.fastq    61144     51990
    ## 144_E0_S18_L001_R1_001.fastq    53714     46734
    ## 144_E4_S19_L001_R1_001.fastq    41686     35201
    ## 144_E8_S20_L001_R1_001.fastq    34947     28743
    ## 144_F0_S21_L001_R1_001.fastq    54554     47488
    ## 144_F4_S22_L001_R1_001.fastq    32800     28274
    ## 144_F8_S23_L001_R1_001.fastq    33312     29691
    ## 144_G0_S24_L001_R1_001.fastq    40935     35533
    ## 144_G4_S25_L001_R1_001.fastq    40109     34624
    ## 144_G8_S26_L001_R1_001.fastq    35610     31334
    ## 144_H0_S27_L001_R1_001.fastq    63711     56558
    ## 144_H4_S28_L001_R1_001.fastq    27892     23924
    ## 144_H8_S29_L001_R1_001.fastq    36860     31648

# Experimental dplyr Extensions for Genomic Range Data

These are some **highly experimental** (e.g. beware, things will change)
functions that make working with genomic range data held in dataframes (or
better, [tibbles](https://github.com/tidyverse/tibble)). If you need more more
elaborate genomic range operations in R, you almost certainly want
Bioconductor's
[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html).


```{R}
dmel_lens <- tibble::tibble(chrom=c("2L", "X", "3L", "4", "2R", "3R"),
                            length = c(23011544L, 22422827L, 24543557L, 
                                       1351857L, 21146708L, 27905053L))

d <- sim_ranges(20e3, chrom_lengths=dmel_lens)
```

`sim_ranges()` automatically appends the chromosome lengths as an attribute to
the dataframe. In general, you'll need to do this on your own:

```{R}
> chrom_lengths(d) <- dmel_lens
> chrom_lengths(d)
# A tibble: 6 × 2
  chrom   length
  <chr>    <int>
1    2L 23011544
2     X 22422827
3    3L 24543557
4     4  1351857
5    2R 21146708
6    3R 27905053
```

We'll add some random data to each element:

```{R}
library(dplyr)
d %>% mutate(data=runif(length(start)))
```

Then, using `create_windows()`, we can add a column of windows to our
dataframe:

```{R}
d %>% create_windows(width=10e3)
```


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

This object `d` is a `gnibble` -- a version of a tibble designed to have the
columns `chrom`, `start`, `end`. (Strand is not implemented yet; this could be
added later, but if you need strand, use GenomicRanges.)

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

Single-table dplyr verbs should work with gnibbles automatically, passing
required metadata in the attributes of this object as needed. `gplyr` overloads
`mutate()`, as this verb does not.  Let's use `mutate()` to add some random
data to each element. 

```{R}
> library(dplyr)
> d <- sim_ranges(20e3, dmel_lens) %>% mutate(x=runif(length(start)))
> class(d)
[1] "gnibble"    "tbl_df"     "tbl"        "data.frame"
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

Then, using `create_windows()`, we can add a column of windows to our
dataframe:

```{R}
> d %>% create_windows(width=10e3)
# A tibble: 20,000 × 5
   chrom start   end       data         window
   <chr> <dbl> <dbl>      <dbl>         <fctr>
1     2L  4370  8415 0.30497986      2L:0-9995
2     2L 19117 28209 0.24789848  2L:9996-19991
3     2L 22177 27374 0.91571345 2L:19992-29987
4     2L 24336 29586 0.40156194 2L:19992-29987
5     2L 31698 40083 0.07644564 2L:29988-39984
6     2L 38233 43566 0.82448720 2L:29988-39984
7     2L 39134 43452 0.98660361 2L:29988-39984
8     2L 60407 60698 0.40600818 2L:59977-69973
9     2L 63615 64942 0.41598095 2L:59977-69973
10    2L 80542 80694 0.12093791 2L:79970-89965
# ... with 19,990 more rows
```
 
With our new `window` column, we can group the data by `chrom`, and `window`,
and then summarize 

```{R}
> library(tidyr)
> d %>% create_windows(1e6) %>% arrange(chrom) %>% 
        group_by(chrom, window) %>% summarize(xmean=mean(x))
Source: local data frame [124 x 3]
Groups: chrom [?]

   chrom             window     xmean
   <chr>             <fctr>     <dbl>
1     2L        2L:0-958813 0.5068166
2     2L  2L:958814-1917627 0.4968698
3     2L 2L:1917628-2876442 0.5111529
4     2L 2L:2876443-3835256 0.4956974
5     2L 2L:3835257-4794070 0.5124214
6     2L 2L:4794071-5752885 0.4825592
7     2L 2L:5752886-6711699 0.5153311
8     2L 2L:6711700-7670514 0.5014180
9     2L 2L:7670515-8629328 0.4714918
10    2L 2L:8629329-9588142 0.5032267
# ... with 114 more rows
>
```



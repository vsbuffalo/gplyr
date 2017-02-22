# Experimental dplyr Extensions for Genomic Range Data

These are some **experimental** (e.g. beware, things will change... including
this uncreative package name) functions that make working with genomic range
data held in dataframes (or better,
[tibbles](https://github.com/tidyverse/tibble)). If you need more more
elaborate genomic range operations in R, you almost certainly want
Bioconductor's
[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html).


```{R}
set.seed(0)
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
   chrom  start    end          x          window
   <chr>  <dbl>  <dbl>      <dbl>          <fctr>
1     2L  12054  14424 0.16426985   2L:9996-19991
2     2L  25353  34091 0.51285709  2L:19992-29987
3     2L  28388  31698 0.73325516  2L:19992-29987
4     2L  34465  44224 0.32336798  2L:29988-39984
5     2L  35375  39163 0.32927339  2L:29988-39984
6     2L  47052  51047 0.58059202  2L:39985-49980
7     2L  72249  73636 0.08094737  2L:69974-79969
8     2L  88559  95516 0.96380628  2L:79970-89965
9     2L  99125 102255 0.31359963  2L:89966-99962
10    2L 101196 108285 0.87439077 2L:99963-109958
# ... with 19,990 more rows
```
 
With our new `window` column, we can group the data by `chrom`, and `window`,
and then summarize 

```{R}
> library(tidyr)
> d %>% create_windows(1e6) %>% 
        group_by(window) %>% 
        summarize_window(xmean=mean(x))
# A tibble: 124 × 3
               window chrom     xmean
               <fctr> <chr>     <dbl>
1         2L:0-958813    2L 0.4964056
2   2L:958814-1917627    2L 0.4932900
3  2L:1917628-2876442    2L 0.5296287
4  2L:2876443-3835256    2L 0.4607771
5  2L:3835257-4794070    2L 0.5052921
6  2L:4794071-5752885    2L 0.5090671
7  2L:5752886-6711699    2L 0.4575886
8  2L:6711700-7670514    2L 0.5106510
9  2L:7670515-8629328    2L 0.4714383
10 2L:8629329-9588142    2L 0.5210464
# ... with 114 more rows
```

`summarize_window()` is like dplyr's `summarize()` except (1) it will group by
`window` if that column is present, and (2) it will carry window variables like
`chrom` forward (and if `wstart` and `wend` are present, those too).

Often, we need to refer to a window's start and end positions. We can extract
these using `separate_window()`. Below, we demonstrate this

```{R}
> dr <- d %>% create_windows(1e6) %>% 
            group_by( window) %>% summarize(xmean=mean(x)) %>%
            separate_window() 
> dr
Source: local data frame [124 x 4]
Groups: chrom [6]

    chrom  wstart    wend     xmean
   <fctr>   <int>   <int>     <dbl>
1      2L       0  958813 0.4964056
2      2L  958814 1917627 0.4932900
3      2L 1917628 2876442 0.5296287
4      2L 2876443 3835256 0.4607771
5      2L 3835257 4794070 0.5052921
6      2L 4794071 5752885 0.5090671
7      2L 5752886 6711699 0.4575886
8      2L 6711700 7670514 0.5106510
9      2L 7670515 8629328 0.4714383
10     2L 8629329 9588142 0.5210464
# ... with 114 more rows
dr <- dr %>% mutate(wcenter=(wstart+wend)/2)
ggplot(dr) + geom_point(aes(x=wcenter, y=xmean)) + facet_wrap(~chrom)
``` 

This workflow is also encapsulated in the convenience function,
`append_wcenter():


```{R}
> d %>% create_windows(1e6) %>% 
        group_by(window) %>% 
        summarize_window(xmean=mean(x)) %>%
        append_wcenter()
# A tibble: 124 × 6
    chrom  wstart    wend             window     xmean   wcenter
   <fctr>   <int>   <int>             <fctr>     <dbl>     <dbl>
1      2L       0  958813        2L:0-958813 0.4964056  479406.5
2      2L  958814 1917627  2L:958814-1917627 0.4932900 1438220.5
3      2L 1917628 2876442 2L:1917628-2876442 0.5296287 2397035.0
4      2L 2876443 3835256 2L:2876443-3835256 0.4607771 3355849.5
5      2L 3835257 4794070 2L:3835257-4794070 0.5052921 4314663.5
6      2L 4794071 5752885 2L:4794071-5752885 0.5090671 5273478.0
7      2L 5752886 6711699 2L:5752886-6711699 0.4575886 6232292.5
8      2L 6711700 7670514 2L:6711700-7670514 0.5106510 7191107.0
9      2L 7670515 8629328 2L:7670515-8629328 0.4714383 8149921.5
10     2L 8629329 9588142 2L:8629329-9588142 0.5210464 9108735.5
# ... with 114 more rows
```

Or, if you want all chromosomes concatenated on a single x-axis, use
`append_wcumpos()`.

```{R}
dr <- d %>% create_windows(1e6) %>% arrange(chrom) %>% 
            group_by(chrom, window) %>% summarize(xmean=mean(x)) %>%
            ungroup() %>%
            separate_window() %>% 
            append_wcumpos()
ggplot(dr) + geom_point(aes(x=wcumpos, y=xmean, color=chrom)) 
``` 



This reprex appears to crash R.
See standard output and standard error for more details.

#### Standard output and error

``` sh
Invalid syntax for chunk options:

load-data echo=FALSE

Please see documentation at https://yihui.org/knitr/options/.

Error in parse(text = code, keep.source = FALSE) : 
  <text>:1:18: unexpected symbol
1: alist( load-data echo
                     ^
```


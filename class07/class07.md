Class 7: R functions and packages
================
Jaya George
4/23/2019

Functions revisited
===================

We will source a file with our functions from last day

``` r
source("http://tinyurl.com/rescale-R")
```

Try out the last day's rescale() function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Try rescale2() function that catches string inputs

``` r
#rescale2(c(1:10, "string"))
```

Finding missing NA values in vectors
====================================

Start with a simple example of a larger problem I am trying to solve

``` r
x <- c(1,2,NA,3,NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

Try putting these together with an AND (TRUE= position with both NA)

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take sum() to find out how many TRUE values we have and thus how many NAs we have in both x and y (because TRUE=1 and FALSE=0)

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

Now I can make this into our first function...

``` r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x, c(NA, 3, NA, 2, NA))
```

    ## [1] 2

Test, test, test (When the length of the vectors aren't the same the first value of the shorter vector is added to the end -&gt; recylcing)

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
```

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
y3 <- c( 1, NA, NA, NA, NA, NA, NA)
both_na(x, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 6

``` r
3==3
```

    ## [1] TRUE

``` r
3 != 2
```

    ## [1] TRUE

``` r
length(x) != length(y2)
```

    ## [1] TRUE

Now lets try both\_na2() function on our different length input vectors

``` r
both_na2 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be the same length")
 }
 sum( is.na(x) & is.na(y) )
}
```

``` r
#both_na2 (x, y2)
```

``` r
which( c(F, F, T, F, T))
```

    ## [1] 3 5

``` r
which(is.na ( c(1,2,NA,4) ))
```

    ## [1] 3

``` r
both_na3 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be vectors of the same length")
 }

 na.in.both <- ( is.na(x) & is.na(y) )
 na.number <- sum(na.in.both)
 na.which <- which(na.in.both)
 message("Found ", na.number, " NA's at position(s):",
 paste(na.which, collapse=", ") )

 return( list(number=na.number, which=na.which) )
}
```

``` r
x <- c(1,2,NA,3,NA)
y <- c(NA, 3, NA, 3, 4)
both_na3(x,y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

Simplify into single vectors

``` r
x <- df1$IDs
y <- df2$IDs

x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
which( x%in%y )
```

    ## [1] 2 3

``` r
x[ x%in%y ]
```

    ## [1] "gene2" "gene3"

``` r
y[ y %in% x]
```

    ## [1] "gene2" "gene3"

A working snippet!

``` r
cbind( x[ x %in% y ], 
       y[ y %in% x ] )
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

Use R studio shortcut 'Code-&gt; Extract function' to turn our snippet into a working function

``` r
gene_intersect <- function(x, y) {
  cbind( x[ x %in% y ], 
         y[ y %in% x ] )
}
```

``` r
gene_intersect(df1$IDs, df2$IDs)
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

Our input $IDs column name may change:add flexibility by allowing the user to specify the gene containing column name

``` r
gene_intersect2 <- function(df1, df2) {
 cbind( df1[ df1$IDs %in% df2$IDs, ],
 df2[ df2$IDs %in% df1$IDs, "exp"] )
}
```

``` r
gene_intersect2(df1, df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

Change column names to something nicer

``` r
gene_intersect3 <- function(df1, df2, gene.colname="IDs") {

 cbind( df1[ df1[,gene.colname] %in%
 df2[,gene.colname], ],
 exp2=df2[ df2[,gene.colname] %in%
 df1[,gene.colname], "exp"] )

}
```

``` r
gene_intersect3(df1, df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

``` r
merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

Grade function -&gt; calculate mean score and drop the lowest single score

``` r
x <- c( 100, 90, 90, 90, 90, 90, 97, 80)
```

``` r
grade <- function(x) {
  (sum(x)- min(x))/ (length(x)-1)
}
```

``` r
grade(x)
```

    ## [1] 92.42857

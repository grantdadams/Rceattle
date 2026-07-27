# Save a ggplot to `file` if a file stem is supplied; always return the plot

Save a ggplot to `file` if a file stem is supplied; always return the
plot

## Usage

``` r
.save_ggplot(p, file = NULL, suffix = "plot", width = 7, height = 6.5)
```

## Arguments

- p:

  A ggplot.

- file:

  File stem (no extension), or `NULL` to skip saving.

- suffix:

  Appended to `file` to form the PNG name.

- width, height:

  Figure size in inches.

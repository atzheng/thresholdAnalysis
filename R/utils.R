## Utility functions
## ===================================================================
## Normalizers
## -------------------------------------------------------------------
#' Does basic checks for edge cases for normalizers
normalizer.check <- function(x){
    if(length(x) <= 1)
        'too small'
    else if(length(unique(x)) <= 1)
        'too few values'
    else
        'valid'
}

#' Centers the vector x around its mean and normalizes by its standard
#' deviation.
standardize <- function(x, na.rm = FALSE){
    if(normalizer.check(x) == 'valid')
        (x - mean(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
    else
        x
}

#' Normalizes the vector x by its maximum value. Note that this
#' doesn't really make sense if all the values of x are negative
max.normalizer <- function(x){
    if(normalizer.check(x) == 'valid')
        x / max(x)
    else
        if(all(x < 0))
            warning('max.normalizer should not be used with vectors of all negative values')
        x
}

#' Subtracts the minimum value of x from x, then normalizes by its
#' range
max.min.normalizer <- function(x){
    if(normalizer.check(x) == 'valid')
        (x - min(x)) / (max(x) - min(x))
    else
        x
}

#' Subtracts the minimum value of x from x, then normalizes by the
#' range from min(x) to its qth quantile
q.min.normalizer <- function(x, q = 0.90){
    if(normalizer.check(x) == 'valid')
        (x - min(x)) / (quantile(x, q) - min(x))
    else
        x
}

# Miscellaneous
# --------------------------------------------------------------------
#' Checks whether a vector x contains only integers. Note that
#' base::is.integer is entirely different.
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)
    abs(x - round(x)) < tol

#' Infers whether a given vector is continuous or discrete based on
#' some basic heuristics: integers are considered discrete, as are
#' non-integer variables where there are less than 20% as many unique
#' values as data points.
is.continuous <- function(x){
    if(all(is.wholenumber(x)) | length(unique(x)) < (length(x) / 5))
        FALSE
    else
        TRUE
}

#' Replaces NAs in x with the most recent non-NA value.
#' Taken from http://stackoverflow.com/questions/7735647/
na.fill.before <- function(x){
    ind = which(!is.na(x))
    if(is.na(x[1]))
        ind = c(1,ind)
    rep(x[ind], times = diff(c(ind, length(x) + 1)))
}

#' Replaces NAs in x with the next non-NA value.
na.fill.after <- function(x)
    rev(na.fill.before(rev(x)))

# Distributions
# --------------------------------------------------------------------
pdf_discrete <- function(x){
    p.x <- unlist(lapply(split(x, x), length)) / length(x)
    return.val <- data.table(x = as(names(p.x), class(x)), p.x = p.x)
    setkey(return.val, x)
    return.val
}

#' Search for an appropriate bandwidth should be implemented here
pdf_continuous <- function(x, bw = 2){
    pdf.x <- density(x, bw = bw)
    return.val <- data.table(x = pdf.x$x, p.x = pdf.x$y)
    setkey(return.val, x)
    return.val
}

d.pdf <- function(x, continuous = NULL, bw = 2){
    if(is.null(continuous))
        continuous <- is.continuous(x)
    if(continuous)
        d.pdf_continuous(x, bw)
    else
        d.pdf_discrete(x)
}

#' A version of the d.pdf function that accepts a precomputed pdf as
#' an argument
d.pdf_precomputed <- function(pdf.x, continuous = NULL, bw = 2){
    if(is.null(continuous))
        continuous <- is.continuous(pdf.x[, x])
    if(continuous)
        d.pdf_continuous(pdf.x, bw)
    else
        d.pdf_discrete(pdf.x)
}

#' Uses the forward finite difference method to compute the first
#' derivative of the pdf of some discrete-valued vector x
d.pdf_discrete <- function(x){
    pdf.x <- pdf_discrete(x)
    return.val <- data.table(
        x = pdf.x[, x],
        d.p.x = c(diff(pdf.x[, p.x]), NA))
    setkey(return.val, x)
    return.val
}

#' For some continuous-valued vector, compute the density using
#' kernel density estimation with gaussian kernels and a manually set
#' bandwidth, and compute the derivative.
#'
#' Method taken from:
#' http://stackoverflow.com/questions/12568715/derivative-of-kernel-density
d.pdf_continuous <- function(x, bw = 2){
    pdf.x <- pdf_continuous(x, bw)
    d.p.x <- vapply(
        unique(x),
        function(y)
            mean(dnorm(x - y, mean = 0, sd = bw) * (x - y)),
        numeric(1))
    return.val <- data.table(x = unique(x), d.p.x = d.p.x)
    setkey(return.val, x)
    return.val
}

create.pdfs <- function(x.list, bw = 2)
    standardize.pdf(lapply(
        x.list,
        function(v){
            if(is.continuous(v))
                pdf_continuous(v, bw)
            else
                pdf_discrete(v)}))

#' Creates a single cdf
#'
#' @param x (numeric) Vector of values from which to generate a cdf.
#' @param inverted (boolean) If TRUE, returns the inverted cdf.
cdf <- function(x, inverted = FALSE){
    x.sorted <- sort(x, decreasing = inverted)
    cum.pct <- seq(1, length(x.sorted))/length(x.sorted)
    return.val <- data.table(x = rev(x.sorted),
                             cum.pct = rev(cum.pct)
                             )[, list(cum.pct = max(cum.pct)), x]
    setkey(return.val, x)
    return.val
}

#' Returns the empirical inverse CDF of any vector x
inv.cdf <- function(x)
    cdf(x, inverted = TRUE)

#' Returns a list of data.tables representing cdfs, where all cdfs are
#' calculated at the same set of values (i.e., standardized).
#'
#' @param x.list (list [numeric]) List of vectors with which to
#' generate cdfs
#' @param inverted (boolean) If TRUE, returns the inverted cdf.
create.cdfs <- function(x.list, inverted = FALSE)
    standardize.cdf(lapply(x.list,
                           function(x) cdf(x, inverted)))

#' Standardizes a list of cdfs such that all are defined at the
#' same values, by filling in previously undefined values with the
#' last non-NA value.
#'
#' Note that since all cdfs should be keyed by their x-values, this
#' function assumes that they are sorted. The case where all cdfs have
#' only one value is assumed to be non-inverted.
#'
#' @param cdf.list A list of data.tables representing cdf's, with
#' columns 'x' and 'cum.pct'
standardize.cdf <- function(cdf.list){
    ## Check whether the cdfs are inverted
    sample.cdf <- Filter(function(cdf) nrow(cdf) > 1, cdf.list)[[1]]
    if(length(sample.cdf) == 0
       | sample.cdf[1, cum.pct] < sample.cdf[2, cum.pct])
        na.fill.fn <- na.fill.before
    else
        na.fill.fn <- na.fill.after

    ## Standardize the cdfs
    all.x <- sort(unlist(Reduce(
        function(init, cdf.x)
            unique(c(init, cdf.x$x)),
        cdf.list,
        init = c())))

    return.val <- lapply(
        cdf.list,
        function(cdf)
            cdf[J(all.x)
                ][, cum.pct := na.fill.fn(cum.pct)
                  ][is.na(cum.pct), cum.pct := 0])
    return.val
}

#' Standardizes a list of pdfs such that all are defined at the
#' same values, by filling in previously undefined values with 0.
#'
#' @param pdf.list A list of data.tables representing pdf's, with
#' columns 'all.x' and 'p.all.x'
standardize.pdf <- function(pdf.list){
    all.x <- sort(unlist(Reduce(function(init, pdf)
        unique(c(init, pdf$x)),
                                pdf.list, init = c())))

    # If continuous, linearly interpolate between points
    if(is.continuous(all.x))
        lapply(pdf.list,
               function(pdf){
                   f <- approxfun(
                       x = pdf$x, y = pdf$p.x, yleft = 0, yright = 0)
                   data.table(x = all.x, p.x = f(all.x))
               })
    else
        lapply(
            pdf.list,
            function(cdf)
                cdf[J(all.x)
                    ][is.na(p.x), p.x := 0])
}


#' ** NOTE: This should be made much more concise with 'standardize'
#'
#' Calculates an inverse cdf for each category, at all unique
#' values of metric.vector.
#'
#' @param category.vector Vector of n-ary categories
#' @param metric.vector Vector for the cdfs
inv.cdf.by.category <- function(category.vector, metric.vector){

    categories <- unique(category.vector)
    metric.values <- sort(unique(metric.vector))

    sorted.order <- order(metric.vector)
    metric.sorted <- metric.vector[sorted.order]
    cat.sorted <- category.vector[sorted.order]

    ## Find the inv.cdf for all categories, at all unique values of
    ## metric.vector
    inv.cdfs <- lapply(
        categories,
        function(category)
            inv.cdf(metric.sorted[cat.sorted == category]))

    inv.cdfs.standardized <- standardize.cdf(inv.cdfs)
    inv.cdfs.standardized$combined <- inv.cdf(metric.vector)

    return.val <- data.table(
        data.frame(lapply(inv.cdfs.standardized,
                          function(cdf) cdf[, cum.pct])))
    names(return.val) <- c(categories, 'combined')
    return.val[, x := metric.values]
    setkey(return.val, x)
    return.val
}

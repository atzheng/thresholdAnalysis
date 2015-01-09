# Alternative Metrics
# ====================================================================
# Jensen-Shannon Divergence
# --------------------------------------------------------------------
jensen.shannon.sqrt <- function(x.A, x.B, continuous = NULL)
    sqrt(jensen.shannon(x.A, x.B, continuous))

jensen.shannon <- function(x.A, x.B, continuous = NULL){
    if(is.null(continuous))
        continuous <- is.continuous(c(x.A, x.B))

    if(continuous)
        jensen.shannon_continuous(x.A, x.B)
    else
        jensen.shannon_discrete(x.A, x.B)
}

jensen.shannon_discrete <- function(x.A, x.B){
    pdf.A <-pdf.discrete(x.A)
    pdf.B <- pdf.discrete(x.B)

    pdfs.std <- standardize.pdf(list(pdf.A, pdf.B))
    pdf.combined <- data.table(
        x = pdfs.std[[1]][, x],
        p.x = (pdfs.std[[1]][, p.x] + pdfs.std[[2]][, p.x]) / 2)

    (kl.divergence_discrete(pdfs.std[[1]], pdf.combined)
     + kl.divergence_discrete(pdfs.std[[2]], pdf.combined)) / 2
}

jensen.shannon_continuous <- function(x.A, x.B){
    cdf.A <- cdf(x.A)
    cdf.B <- cdf(x.B)

    cdfs.std <- standardize.cdf(list(cdf.A, cdf.B))
    cdf.combined <- data.table(
        x = cdfs.std[[1]][, x],
        cum.pct = (cdfs.std[[1]][, cum.pct]
                   + cdfs.std[[2]][, cum.pct]) / 2)

    (kl.divergence_continuous(cdf.A, cdf.combined)
     + kl.divergence_continuous(cdf.B, cdf.combined)) / 2
}

#' Computes the Kullback-Leibler divergence of two discrete
#' distributions, given their pdfs.
#'
#' Still need to implement smoothing instead of na.rm
kl.divergence_discrete <- function(pdf.A, pdf.B, smoothing = 1)
    sum(pdf.A[, p.x] * log(pdf.A[, p.x] / pdf.B[, p.x]), na.rm = TRUE)

#' Computes the Kullback-Leibler divergence of two continuous
#' distributions, given their empirical cdfs. This uses the
#' approximation specified in "Kullback-Leibler Divergence Estimation
#' of Continuous Distributions" (Perez-Cruz, 2008)
kl.divergence_continuous <- function(cdf.A, cdf.B){
    fun.A <- approxfun(cdf.A[, x], cdf.A[, cum.pct])
    fun.B <- approxfun(cdf.B[, x], cdf.B[, cum.pct])
    all.x <- sort(unique(c(cdf.A$x, cdf.B$x)))
    epsilon <- min(diff(all.x)) / 2

    p.q <- ((fun.A(all.x[-1]) - fun.A(all.x[-1] - epsilon))
            / (fun.B(all.x[-1]) - fun.B(all.x[-1] - epsilon)))
    mean(log(p.q), na.rm = TRUE)
}

# CDF Separation
# --------------------------------------------------------------------
separation.cdf <- function(x.A, x.B){
    cdfs <- create.cdfs(list(x.A, x.B), inverted = TRUE)
    cdf.A <- cdfs[[1]]
    cdf.B <- cdfs[[2]]
    x <- cdf.A[, x]

    x.distance <- diff(x)
    separation <- sum(abs(cdf.A[-1, cum.pct]
                          - cdf.B[-1, cum.pct]) * x.distance)
    return(separation)
}

#' Returns a metric of separation between two functions, defined as
#' the sum of the area bounded on both sides by the two functions.
#'
#' @param x (numeric) Vector of values for which f(x) and g(x) are
#' computed
#' @param f.x (numeric) Vector of f(x) values at all values of x
#' @param g.x (numeric) Vector of g(x) values at all values of x
#'
#' @return (numeric.scalar) The value of the separation metric
separation.metric <- function(x, f.x, g.x){
    x.distance <- diff(x)
    separation <- sum(abs(g.x[-1] - f.x[-1]) * x.distance)
    return(separation)
}

# Linearly interpolated pdf separation (Inefficient Implementation)
# --------------------------------------------------------------------
anomaly.metric <- function(x.A, x.B, bin.size = NULL, h = 1){
    cdfs <- create.cdfs(list(x.A, x.B), inverted = TRUE)
    cdf.A <- cdfs[[1]]
    cdf.B <- cdfs[[2]]
    anomaly.metric_helper(cdf.A[, x],
                          cdf.A[, cum.pct],
                          cdf.B[, cum.pct],
                          bin.size,
                          h)
}

#' Returns a metric of anomaly detection between two functions,
#' defined as the sum of the area bounded on both sides by the
#' derivatives of the two functions.
#'
#' @param x (numeric) Vector of values for which f(x) and g(x) are
#' computed
#' @param f.x (numeric) Vector of f(x) values at all values of x
#' @param g.x (numeric) Vector of g(x) values at all values of x
#' @param bin.size (numeric.scalar) Intervals at which to compute the
#' first derivative of f.x and g.x
#' @param h (numeric.scalar) Peturbation parameter to the finite
#' differences method
#'
#' @return (numeric.scalar) The value of the anomaly metric
anomaly.metric_helper <- function(x, f.x, g.x, bin.size = NULL, h = 1){
    d.dx <- lapply(list(f.x, g.x),
                   function(h.x)
                       anomaly.curve(x, h.x, bin.size, h))
    x.bins <- d.dx[[1]][, x]
    separation.metric(x.bins, d.dx[[1]][, df.dx], d.dx[[2]][, df.dx])
}

#' Generates the first derivative curve of the step function f.x at
#' values of x defined by
#'
#' @param x (numeric) Vector of values for which f(x) is computed
#' @param f.x (numeric) Vector of f(x) values at all values of x
anomaly.curve <- function(x, f.x, bin.size = NULL, h = 1){
    if(is.null(bin.size)){
        bin.size <- max(1, max(x) / 200)
    }
    x.bins <- seq(min(x), max(x), bin.size)
    df.dx <- unlist(lapply(
        x.bins,
        function(thresh)
            finite.differences(x, f.x, thresh, 1, h)))
    data.table(x = x.bins, df.dx)
}

# Finite differences method
# --------------------------------------------------------------------
#' Uses the central finite difference method to calculate the
#' second derivative of stability.vector's empirical inverse cdf at
#' threshold. Note that if h is not large enough (required size
#' varies inversely with stability.vector's size), stability will
#' be Inf for most values.
#'
#' stability = \frac{1}{\frac{d (cdf(depth, x))}{dx} + 1}
stability.metric <- function(x, f.x, threshold, h = 1){
    return(finite.differences(x, f.x,  threshold, 2, h))
}

#' Calculates the <level> derivative of metric.vector's cdf using the
#' central finite differences method.
#'
#' @param x values at which f.x is computed
#' @param f.x values of f(x) computed at all values of x
#' @param threshold value at which to compute the derivative
#' @param level derivative level (1 -> first derivative, 2 -> second
#' derivative)
#' @param h Perturbation over which to compute the derivative
#'
#' TODO: This doesn't handle the case where there is no x between
#  thresh and thresh + h
finite.differences <- function(x,
                               f.x,
                               threshold,
                               level = 2,
                               h = 1){
    min.bound <- min(x)
    max.bound <- max(x)

    lower.bound <- threshold - h
    upper.bound <- threshold + h
    normalizing.factor <- 2 * h

    ## -- Handle edge cases:
    ## (1) If the threshold is outside of the range, return NA
    if(threshold < min.bound
       | threshold > max.bound
       | (lower.bound < min.bound & upper.bound > max.bound)){
        return(NA)
    }
    ## (2) If the threshold is within the range but the lower or upper
    ##     bound is not, then use the forward or backward finite
    ##     difference method.
    if(lower.bound < min.bound)
        lower.bound <- threshold
    if(upper.bound > max.bound)
        upper.bound <- threshold
    ## -- Calculate the derivative
    ## (1) If we only need the first derivative, calculate and return
    if(level == 1){
        ## Since f(x) is always a step function, f(x') for any x'
        ## not in x is f(x''), where x'' = max(x) s.t. x <= x'
        min.x <- max(x[x <= lower.bound])
        max.x <- min(x[x >= upper.bound])
        return((min(f.x[x == max.x]) - max(f.x[x == min.x]))
               / (max.x - min.x))
    }
    ## (2) If we need the nth derivative, find it recursively
    else{
        necessary.x <- unique(
            pmin(max.bound,
                 pmax(min.bound,
                      c(seq(0, level), -seq(level)) * h + threshold)))
        df.dx <- unlist(lapply(
            necessary.x,
            function(thresh)
                finite.differences(
                    x, f.x, thresh, 1, h)))
        return(
            finite.differences(
                necessary.x, df.dx, threshold, level - 1, h))
    }
}


# Separation
# --------------------------------------------------------------------
#' Computes the absolute y-axis distance between two probability
#' density functions. In the continuous case, these pdf's are
#' approximated using kernel density estimation with some default
#' parameters; these can be altered by providing an alternate density
#' estimation function in the function parameters.
#'
#' @param x.A (numeric) Vector of values for computing the first pdf
#' @param x.B (numeric) Vector of valeus for computing the second pdf
#' @param continuous (boolean scalar) Specifies whether the data is
#' continuous or discrete; if NULL (default), the function will infer
#' this from the data.
#' @param density.fn (function) Function to use to compute the density
#' of a continuous variable (only necessary if continuous). This must
#' be of the form f(x) -> list(x = x, y = p.x)
#' @param subdivisions (integer) Only necessary if continuous; number
#' of subdivisions to use to compute the integral of a continuous pdf.
separation.pdf <- function(x.A,
                           x.B,
                           continuous = NULL,
                           density.fn = density,
                           subdivisions = max(length(x.A),
                                              length(x.B)) * 2){
    # Infer whether the data is continuous or discrete if not
    # specified
    if(is.null(continuous))
        continuous <- is.continuous(c(x.A, x.B))

    # Discrete case
    if(continuous)
        separation.pdf_continuous(x.A, x.B, density.fn, subdivisions)
    else
        separation.pdf_discrete(x.A, x.B)
}

#' Helper function for pdf.separation; handles the discrete case
separation.pdf_discrete <- function(x.A, x.B){
    # Compute the pdfs of x.A and x.B
    pdfs <- create.pdfs(list(x.A, x.B))
    pdf.A <- pdfs[[1]]
    pdf.B <- pdfs[[2]]

    sum(abs(pdf.A[, p.x] - pdf.B[, p.x]))
}

#' Helper function for pdf.separation; handles the continuous case.
separation.pdf_continuous <- function(
    x.A, x.B, density.fn, subdivisions){

    pdfs <- lapply(
        list(x.A, x.B),
        function(v){
            ## Calculate initial pdf
            pdf.x <- density.fn(v)

            ## Truncate pdf to known range and renormalize
            in.range <- pdf.x$x >= min(v) & pdf.x$x <= max(v)
            pdf.x.norm <- pdf.x$x[in.range]
            pdf.y <- pdf.x$y[in.range]
            approx.pdf <- approxfun(
                pdf.x.norm, pdf.y, yleft = 0, yright = 0)
            normalizer <- integrate(
                approx.pdf, min(v), max(v),
                subdivisions = subdivisions)$value
            pdf.y.norm <- pdf.y / normalizer
            approx.pdf.norm <- approxfun(
                pdf.x.norm, pdf.y.norm, yleft = 0, yright = 0)
        })

    integrate(function(x) abs(pdfs[[1]](x) - pdfs[[2]](x)),
              min(c(x.A, x.B)),
              max(c(x.A, x.B)),
              subdivisions = subdivisions)
}

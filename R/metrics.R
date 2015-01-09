# Thresholding Metrics
# ====================================================================
# Note that many of these metrics share several of the computations
# (i.e. calculation of the pdf), so it's possible to implement these
# much more efficiently; however, for the sake of clarity and
# modularity, we separate the computations here.

# Stability
# --------------------------------------------------------------------
#' Computes the first derivative of the pdf of x at all values of x
stability <- function(x, thresholds, continuous = NULL, bw = 2)
    d.pdf(x, continuous, bw)[thresholds, d.p.x]

## Separation
## -------------------------------------------------------------------
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
                           subdivisions = 20000){
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
                subdivisions = subdivisions,
                abs.tol = 0.01)$value
            pdf.y.norm <- pdf.y / normalizer
            approx.pdf.norm <- approxfun(
                pdf.x.norm, pdf.y.norm, yleft = 0, yright = 0)
        })

    #' When bootstrapping, the integral occasionally fails with error
    #' "the integral is probably divergent"; unsure what's causing this,
    #' so just return 0 for now
    tryCatch(
        integrate(function(x) abs(pdfs[[1]](x) - pdfs[[2]](x)),
                  min(c(x.A, x.B)),
                  max(c(x.A, x.B)),
                  subdivisions = subdivisions,
                  abs.tol = 0.01)$value,
        error = function(e) 0)
}

## Anomaly
## -------------------------------------------------------------------
anomaly.pdf <- function(x.A, x.B, continuous = NULL){
    if(is.null(continuous))
        continuous <- is.continuous(c(x.A, x.B))
    if(continuous)
        anomaly.pdf_continuous(x.A, x.B)
    else
        anomaly.pdf_discrete(x.A, x.B)
}

anomaly.pdf_discrete <- function(x.A, x.B){
    pdfs <- create.pdfs(list(x.A, x.B))
    d.pdf.A <- diff(pdfs[[1]][, p.x])
    d.pdf.B <- diff(pdfs[[2]][, p.x])
    sum(abs(d.pdf.A - d.pdf.B))
}

anomaly.pdf_continuous <- function(x.A, x.B){
    d.pdf.A <- d.pdf_continuous(x.A)
    d.pdf.B <- d.pdf_continuous(x.B)
    f.A <- approxfun(d.pdf.A[, x], d.pdf.A[, d.p.x])
    f.B <- approxfun(d.pdf.B[, x], d.pdf.B[, d.p.x])

    x <- unique(c(x.A, x.B))
    N <- length(x)
    d.x <- diff(x)
    sum(d.x * abs(f.A(x[-N]) - f.B(x[-N])), na.rm = TRUE)
}

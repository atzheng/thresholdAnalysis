#' This is the general purpose function for the thresholding
#' algorithm.
#'
#' @param data A data.table or data.frame on which to apply the
#' algorithm
#' @param formula A formula specifying the variables in the
#' algorithm. This must be of the form <categorizing variable> ~
#' <noncategorizing variables>
#' @param min.cluster.size The minimum number of data points in any
#' cluster, as a proportion of the total number of data points.
threshold.analysis <- function(data,
                               formula,
                               min.cluster.size = 0.2,
                               num.thresholds = 1,
                               parallel = FALSE){
    if(!is.data.table(data))
        data <- as.data.table(data)

    cat.var <- all.vars(formula)[1]
    data.cat.var <- data[, cat.var, with = FALSE][[1]]
    cat.cdf <- inv.cdf(data.cat.var)
    candidate.thresholds <- (
        cat.cdf[cum.pct >= min.cluster.size
                & cum.pct <= (1 - min.cluster.size), x])

    metrics <- thresholding.metrics(
        data,
        formula,
        thresholds = candidate.thresholds,
        cat.metric.fns = list(stability = stability),
        noncat.metric.fns = list(
            separation = function(x.A, x.B)
                separation.pdf(x.A, x.B, density.fn = function(x)
                    density(x, bw = 3)),
            anomaly = anomaly.pdf),
        parallel = parallel)

    weighted.thresholds(
        metrics[, threshold],
        metrics[, stability],
        metrics[, -which(names(metrics)
                %in% c('threshold', 'stability')),
                with = FALSE],
        num.thresholds = num.thresholds)
}

#' Calculates metrics on the data for a range of thresholds.
#'
#' @param data The dataset on which to calculate metrics
#' @param formula The formula for the metrics, specified as
#' <categorizing.variable> ~ <non.categorizing.variables>
#' @param thresholds A vector of thresholds
#' @param cat.metric.fns A named list of metrics (functions) to
#' calculate over the categorizing variable. These must have the form:
#' f(categorizing.variable.vector, threshold) -> float
#' @param noncat.metric.fns A named list of metrics (functions) to
#' calculate over the noncategorizing variables. Assuming we have
#' categories A and B, these functions must have the form:
#' f(x, cdf.A(x), cdf.B(x)) -> float
#' @param parallel Boolean indicating whether or not to run
#' calculations in parallel by replacing lapply with mclapply
thresholding.metrics <- function(data,
                                 formula,
                                 thresholds,
                                 cat.metric.fns = list(
                                     stability = stability),
                                 noncat.metric.fns = list(
                                     separation = separation.pdf,
                                     anomaly = anomaly.metric),
                                 parallel = FALSE){

    if(parallel) apply.fn <- mclapply
    else apply.fn <- lapply

    ## Extract the variable names from the given formula
    var.names <- all.vars(formula)
    var.list <- apply.fn(var.names,
                         function(var)
                             data[, eval(parse(text = var))])
    names(var.list) <- var.names

    results <- apply.fn(
        thresholds,
        function(thresh)
            thresholding.metrics_helper(var.list[[1]],
                                        var.list[-1],
                                        thresh,
                                        cat.metric.fns,
                                        noncat.metric.fns))

    return.val <- Reduce(rbind, results, init = data.table())
    return.val[, threshold := thresholds]
    return(return.val)
}

#' Helper function for thresholding.metrics that calculates metrics
#' for a single threshold. Must be provided with the variables
#' in vector form.
thresholding.metrics_helper <- function(cat.var,
                                        noncat.vars,
                                        threshold,
                                        cat.metric.fns,
                                        noncat.metric.fns){
    cat.results <- lapply(cat.metric.fns,
                          function(f)
                              f(cat.var,
                                threshold))
    names(cat.results) <- names(cat.metric.fns)

    noncat.results <- unlist(lapply(
        noncat.metric.fns,
        function(metric.fn)
            lapply(noncat.vars,
                   function(var)
                       metric.fn(var[cat.var <= threshold],
                                 var[cat.var > threshold]))))

    return(c(cat.results, noncat.results))
}

# Optimization: Weighted Thresholds
# ====================================================================
#' Calculates the optimal thresholds over a range of weights
#' (lambdas).
weighted.thresholds <- function(thresholds,
                                primary.objective,
                                secondary.objectives,
                                normalizer = max.min.normalizer,
                                num.thresholds = 1){

    ## The first component of the objective is just the primary
    ## objective; we need to calculate the second component
    secondary.normalizer <- length(secondary.objectives)
    secondary.mtx <- do.call(cbind, secondary.objectives)
    if(!is.null(normalizer)){
        secondary.mtx <- apply(secondary.mtx, 2, normalizer)
        primary.objective <- normalizer(primary.objective)
    }
    secondary.objective <- apply(secondary.mtx, 1, function(v)
        sum(v)/secondary.normalizer)

    ## Calculate the best threshold for all lambdas
    optimal.thresholds <- optimize.over.lambdas(
        thresholds, primary.objective, secondary.objective)
    head(optimal.thresholds[order(
        lambda_range, decreasing = TRUE), threshold],  num.thresholds)
}

#' Calculates the optimal thresholds over lambdas from 0 to 1, as well
#' as the size of the range of lambdas for which they are optimal.
#'
#' The optimizer is implemented in Python. For now, this function
#' writes the data to file as a kludge because rPython isn't working
#' properly; in the near future, it should be possible to use
#' rPython to interface directly with the Python code instead.
optimize.over.lambdas <- function(thresholds, f1, f2){
    expand.fname <- function(fname){
        system.file('src', fname, package = 'ThresholdAnalysis')
    }
    write.csv(data.frame(t = thresholds, f1 = f1, f2 = f2),
              file = expand.fname('tmp_in.csv'),
              row.names = FALSE)
    syscmd <- paste("python2.7",
                    expand.fname("lp_optimizer.py"),
                    expand.fname("tmp_in.csv"),
                    expand.fname("tmp_out.csv"))
    system(syscmd)
    return_val <- fread(expand.fname('tmp_out.csv'))
    ## system(paste("rm",
    ##              expand.fname("tmp_in.csv"),
    ##              expand.fname("tmp_out.csv")))
    return_val
}

bootstrap.thresholder <- function(data, index)
    nanohub.thresholder(data[index], num.thresholds = 1)

nanohub.thresholder <- function(data, num.thresholds){
    depth.cdf <- inv.cdf(data[, depth])
    candidate.thresholds <- (
        depth.cdf[cum.pct >= 0.2 & cum.pct <= 0.8, x])

    metrics <- thresholding.metrics(
        data,
        depth.trunc ~ breadth.trunc + length.max.trunc,
        thresholds = candidate.thresholds,
        cat.metric.fns = list(stability = stability),
        noncat.metric.fns = list(
            separation = function(x.A, x.B)
                separation.pdf(x.A, x.B, density.fn = function(x)
                    density(x, bw = 3)),
            anomaly = anomaly.pdf))

    weighted.thresholds(metrics[, threshold],
                        metrics[, stability],
                        metrics[, -which(names(metrics)
                                         %in% c('threshold', 'stability')),
                                with = FALSE],
                        num.thresholds = num.thresholds)
}

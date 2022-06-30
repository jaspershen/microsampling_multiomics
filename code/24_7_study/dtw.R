x = 1:10
y = 1:20

library(dtw)

dtw <-
  function (x,
            y = NULL,
            dist.method = "Euclidean",
            step.pattern = symmetric2,
            window.type = "none",
            keep.internals = FALSE,
            distance.only = FALSE,
            open.end = FALSE,
            open.begin = FALSE,
            ...)
  {
    lm <- NULL
    if (is.null(y)) {
      if (!missing(dist.method))
        stop("Argument dist.method does not make sense with a local cost matrix")
      if (!is.matrix(x))
        stop("Single argument requires a pre-computed local cost matrix")
      lm <- x
    }else {
      x <- as.matrix(x)
      y <- as.matrix(y)
      if ((ncol(x) == 1 || ncol(y) == 1) && !missing(dist.method)){
        warning(
          "Argument dist.method does not usually make a difference with single-variate timeseries"
        )
      }
        
      if (!proxy::pr_DB$entry_exists(dist.method)){
        stop("dist.method should be one of the method names supported by proxy::dist()")
      }
        
      lm <- proxy::dist(x, y, method = dist.method)
    }
    wfun <- .canonicalizeWindowFunction(window.type)
    dir <- step.pattern
    norm <- attr(dir, "norm")
    if (!is.null(list(...)$partial)) {
      warning("Argument `partial' is obsolete. Use `open.end' instead")
      open.end <- TRUE
    }
    n <- nrow(lm)
    m <- ncol(lm)
    if (open.begin) {
      if (is.na(norm) || norm != "N") {
        stop(
          "Open-begin requires step patterns with 'N' normalization (e.g. asymmetric, or R-J types (c)). See papers in citation()."
        )
      }
      lm <- rbind(0, lm)
      np <- n + 1
      precm <- matrix(NA, nrow = np, ncol = m)
      precm[1,] <- 0
    } else {
      precm <- NULL
      np <- n
    }
    gcm <-
      globalCostMatrix(
        lm,
        step.matrix = dir,
        window.function = wfun,
        seed = precm
        # ...
      )
    gcm$N <- n
    gcm$M <- m
    gcm$call <- match.call()
    gcm$openEnd <- open.end
    gcm$openBegin <- open.begin
    gcm$windowFunction <- wfun
    lastcol <- gcm$costMatrix[np,]
    if (is.na(norm)) {
      
    }
    else if (norm == "N+M") {
      lastcol <- lastcol / (n + (1:m))
    }
    else if (norm == "N") {
      lastcol <- lastcol / n
    }
    else if (norm == "M") {
      lastcol <- lastcol / (1:m)
    }
    gcm$jmin <- m
    if (open.end) {
      if (is.na(norm)) {
        stop("Open-end alignments require normalizable step patterns")
      }
      gcm$jmin <- which.min(lastcol)
    }
    gcm$distance <- gcm$costMatrix[np, gcm$jmin]
    if (is.na(gcm$distance)) {
      stop("No warping path exists that is allowed by costraints")
    }
    if (!is.na(norm)) {
      gcm$normalizedDistance <- lastcol[gcm$jmin]
    }
    else {
      gcm$normalizedDistance <- NA
    }
    if (!distance.only) {
      mapping <- backtrack(gcm)
      gcm <- c(gcm, mapping)
    }
    if (open.begin) {
      gcm$index1 <- gcm$index1[-1] - 1
      gcm$index1s <- gcm$index1s[-1] - 1
      gcm$index2 <- gcm$index2[-1]
      gcm$index2s <- gcm$index2s[-1]
      lm <- lm[-1,]
      gcm$costMatrix <- gcm$costMatrix[-1,]
      gcm$directionMatrix <- gcm$directionMatrix[-1,]
    }
    if (!keep.internals) {
      gcm$costMatrix <- NULL
      gcm$directionMatrix <- NULL
    }
    else {
      gcm$localCostMatrix <- lm
      if (!is.null(y)) {
        gcm$query <- x
        gcm$reference <- y
      }
    }
    class(gcm) <- "dtw"
    return(gcm)
  }
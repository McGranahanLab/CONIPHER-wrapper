### Copy of functions from sciClone R package version 1.1.0 added by Ariana Huebner
### To cite the pacakge please use: Chris Miller, Brian White and Nathan Dees (2016). sciClone: Detect subclones from genomic sequence data. R package version 1.1.0.

initScClass <- function() {
  setClass("scObject", representation(clust = "list", densities = "list", dimensions = "numeric",
                                      marginalClust = "list", sampleNames = "character", vafs.1d = "list",
                                      vafs.merged = "data.frame", ellipse.metadata = "list", purities = "numeric"))
}

initScClass()

sciClone <- function (vafs, copyNumberCalls = NULL, regionsToExclude = NULL, 
    sampleNames, minimumDepth = 100, clusterMethod = "bmm", clusterParams = "no.apply.overlapping.std.dev.condition", 
    cnCallsAreLog2 = FALSE, useSexChrs = TRUE, doClustering = TRUE, 
    verbose = TRUE, copyNumberMargins = 0.25, maximumClusters = 10, 
    annotation = NULL, doClusteringAlongMargins = TRUE, plotIntermediateResults = 0) {
    if (verbose) {
        print("checking input data...")
    }
    dimensions = NULL
    if (is.data.frame(vafs)) {
        dimensions = 1
        vafs = list(vafs)
        copyNumberCalls = list(copyNumberCalls)
    }
    else if (is.list(vafs)) {
        dimensions = length(vafs)
    }
    else {
        stop("input param vafs must be either a data frame (for 1-sample clustering), or a list of data frames (for multi-sample clustering)")
    }
    if (missing(sampleNames)) {
        stop("sampleNames is a required parameter")
    }
    if (dimensions > 1) {
        if (length(sampleNames) != dimensions) {
            stop(paste("the number of sample names (", length(sampleNames), 
                ") does not equal the number of input samples (", 
                dimensions, ")", sep = ""))
        }
    }
    if (!(is.null(copyNumberCalls))) {
        if (length(copyNumberCalls) != dimensions) {
            stop(paste("the number of input copy number calls(", 
                length(copyNumberCalls), ") does not equal the number of input samples (", 
                dimensions, ")", sep = ""))
        }
    }
    else {
        print("No copy number files specified. Assuming all variants have a CN of 2.")
        copyNumberCalls = vector("list", dimensions)
    }
    doPurityScaling = FALSE
    purities = NULL
    if (!(is.null(purities))) {
        if (length(purities) != dimensions) {
            stop(paste("the number of input purities calls(", 
                length(purities), ") does not equal the number of input samples (", 
                dimensions, ")\nEither provide a purity for each sample, or set purities to NULL and it will be estimated for you", 
                sep = ""))
        }
    }
    if (!is.null(regionsToExclude)) {
        if (is.data.frame(regionsToExclude)) {
            regionsToExclude = list(regionsToExclude)
        }
    }
    densityData = NULL
    for (i in 1:dimensions) {
        vafs[[i]] = cleanAndAddCN(vafs[[i]], copyNumberCalls[[i]], 
            i, cnCallsAreLog2, regionsToExclude, useSexChrs, 
            minimumDepth, copyNumberMargins)
        if (is.null(densityData)) {
            densityData = list(getDensity(vafs[[i]], copyNumberMargins, 
                minimumDepth))
        }
        else {
            densityData = c(densityData, list(getDensity(vafs[[i]], 
                copyNumberMargins, minimumDepth)))
        }
        if (is.null(densityData[[i]]$densities[[2]])) {
            cat(paste("can't do clustering - no copy number 2 regions to operate on in sample", 
                i, "\n"))
            if (doClustering == TRUE) {
                return(NULL)
            }
        }
    }
    vafs.merged = vafs[[1]]
    if (dimensions > 1) {
        for (i in 2:dimensions) {
            vafs.merged = merge(vafs.merged, vafs[[i]], by.x = c(1, 
                2), by.y = c(1, 2), suffixes = c(i - 1, i), all.x = TRUE, 
                all.y = TRUE)
        }
    }
    refcols = grep("^ref", names(vafs.merged))
    varcols = grep("^var", names(vafs.merged))
    vafcols = grep("^vaf", names(vafs.merged))
    depthcols = grep("^depth", names(vafs.merged))
    cncols = grep("^cn", names(vafs.merged))
    cleancncols = grep("^cleancn", names(vafs.merged))
    for (i in c(vafcols, refcols, varcols, depthcols)) {
        vafs.merged[is.na(vafs.merged[, i]), i] = 0
    }
    names(vafs.merged)[refcols] = paste(sampleNames, ".ref", 
        sep = "")
    names(vafs.merged)[varcols] = paste(sampleNames, ".var", 
        sep = "")
    names(vafs.merged)[vafcols] = paste(sampleNames, ".vaf", 
        sep = "")
    names(vafs.merged)[depthcols] = paste(sampleNames, ".depth", 
        sep = "")
    names(vafs.merged)[cncols] = paste(sampleNames, ".cn", sep = "")
    names(vafs.merged)[cleancncols] = paste(sampleNames, ".cleancn", 
        sep = "")
    adequateDepth = rep(1, length(vafs.merged[, 1]))
    for (i in 1:length(vafs.merged[, 1])) {
        if (length(which(vafs.merged[i, depthcols] >= minimumDepth)) < 
            length(depthcols)) {
            adequateDepth[i] = 0
        }
    }
    vafs.merged$adequateDepth = adequateDepth
    cnNeutral = rep(1, length(vafs.merged[, 1]))
    for (i in 1:length(vafs.merged[, 1])) {
        for (j in cleancncols) {
            if (is.na(vafs.merged[i, j])) {
                cnNeutral[i] = 0
            }
            else if (vafs.merged[i, j] != 2) {
                cnNeutral[i] = 0
            }
        }
    }
    vafs.merged.cn2 = vafs.merged[(as.logical(cnNeutral) & as.logical(adequateDepth)), 
        ]
    if (length(vafs.merged.cn2[, 1]) < 1) {
        print("ERROR: no sites are copy number neutral and have adequate depth in all samples")
        return(NULL)
    }
    cat(paste(length(vafs.merged.cn2[, 1]), "sites (of", length(vafs.merged[, 
        1]), "original sites) are copy number neutral and have adequate depth in all samples\n"))
    cat(paste(nrow(vafs.merged[!as.logical(cnNeutral), ])), "sites (of", 
        nrow(vafs.merged), "original sites) were removed because of copy-number alterations\n")
    cat(paste(nrow(vafs.merged[!as.logical(adequateDepth), ])), 
        "sites (of", nrow(vafs.merged), "original sites) were removed because of inadequate depth\n")
    cat(paste(nrow(vafs.merged[!as.logical(cnNeutral) | !as.logical(adequateDepth), 
        ])), "sites (of", nrow(vafs.merged), "original sites) were removed because of copy-number alterations or inadequate depth\n")
    nonvafcols <- (1:length(names(vafs.merged)))[!((1:length(names(vafs.merged))) %in% 
        vafcols)]
    vafs.merged.orig <- vafs.merged
    vafs.merged.cn2.orig <- vafs.merged.cn2
    vafs.matrix = as.matrix(vafs.merged.cn2[, vafcols])
    vars.matrix = as.matrix(vafs.merged.cn2[, varcols])
    refs.matrix = as.matrix(vafs.merged.cn2[, refcols])
    vafs.matrix = vafs.matrix/100
    if (is.null(purities)) {
        if (doPurityScaling) {
            purities = c()
            print("estimating purity")
            clust = clusterVafs(vafs.merged.cn2, vafs.matrix, 
                vars.matrix, refs.matrix, maximumClusters, clusterMethod, 
                clusterParams, samples = length(purities), plotIntermediateResults = 0, 
                verbose = 0)
            if (is.null(clust[[1]])) {
                print("WARNING: couldn't estimate and correct for purity, will cluster using input vafs")
            }
            else {
                for (i in 1:dimensions) {
                  purities[i] = round(max(sort(t(clust[["cluster.means"]])[, 
                    i])) * 2 * 100, 2)
                  vafs.tmp = cbind(vafs.merged.cn2, cluster = clust$cluster.assignments)
                  vafcols = grep("vaf$", names(vafs.tmp))
                  print(paste("purity of sample", sampleNames[i], 
                    "is estimated to be", purities[[i]]))
                  vafcols = grep("vaf$", names(vafs.merged))
                  vafs.merged[, vafcols[i]] = (vafs.merged[, 
                    vafcols[i]]/purities[i]) * 100
                  vafcols = grep("vaf$", names(vafs.merged.cn2))
                  vafs.merged.cn2[, vafcols[i]] = (vafs.merged.cn2[, 
                    vafcols[i]]/purities[i]) * 100
                  vafs.matrix[, i] = (vafs.matrix[, i]/purities[i]) * 
                    100
                }
            }
        }
        else {
            purities = rep(100, dimensions)
        }
    }
    else {
        if (doPurityScaling) {
            for (i in 1:dimensions) {
                vafs.merged[, vafcols[i]] = (vafs.merged[, vafcols[i]]/purities[i]) * 
                  100
                vafs.merged.cn2[, vafcols[i]] = (vafs.merged.cn2[, 
                  vafcols[i]]/purities[i]) * 100
                vafs.matrix[, i] = (vafs.matrix[, i]/purities[i]) * 
                  100
            }
        }
    }
    marginalClust = list()
    vafs.1d = list()
    if (dimensions == 1) {
        doClusteringAlongMargins <- FALSE
    }
    if (doClustering == FALSE) {
        doClusteringAlongMargins <- FALSE
    }
    if (doClusteringAlongMargins == TRUE) {
        print("clustering each dimension independently")
        for (i in 1:dimensions) {
            marginalClust[[i]] = clusterVafs(NULL, vafs.matrix[, 
                i, drop = FALSE], vars.matrix[, i, drop = FALSE], 
                refs.matrix[, i, drop = FALSE], maximumClusters, 
                clusterMethod, clusterParams, FALSE)
            if (verbose) {
                print(paste("finished 1d clustering", sampleNames[i], 
                  "..."))
            }
            numClusters = max(marginalClust[[i]]$cluster.assignments, 
                na.rm = T)
            print(paste("found", numClusters, "clusters using", 
                clusterMethod, "in dimension", sampleNames[i]))
            print(marginalClust[[i]]$cluster.means)
            vafs.1d.merged.cn2 = cbind(vafs.merged.cn2.orig, 
                cluster = marginalClust[[i]]$cluster.assignments)
            vafs.1d.merged = merge(vafs.merged.orig, vafs.1d.merged.cn2, 
                by.x = c(1:length(vafs.merged.orig)), by.y = c(1:length(vafs.merged.orig)), 
                all.x = TRUE)
            vafs.1d.merged = vafs.1d.merged[order(vafs.1d.merged[, 
                1, drop = FALSE], vafs.1d.merged[, 2, drop = FALSE]), 
                ]
            vafs.1d[[i]] = vafs.1d.merged
        }
    }
    clust = list(NULL)
    if (doClustering) {
        if (verbose) {
            print("clustering...")
        }
        clust = clusterVafs(vafs.merged.cn2, vafs.matrix, vars.matrix, 
            refs.matrix, maximumClusters, clusterMethod, clusterParams, 
            samples = length(purities), plotIntermediateResults = plotIntermediateResults, 
            verbose = 0)
        if (is.null(clust[[1]])) {
            print("Warning: no clusters, returning NULL")
            return(NULL)
        }
        if (verbose) {
            print("finished clustering full-dimensional data...")
        }
        clust = reorderClust(clust)
    }
    numClusters = 0
    if (!(is.null(clust[[1]]))) {
        numClusters = max(clust$cluster.assignments, na.rm = T)
        vafs.merged.cn2 = cbind(vafs.merged.cn2, cluster = clust$cluster.assignments)
        vafs.merged.cn2 = cbind(vafs.merged.cn2, cluster.prob = clust$cluster.probabilities)
        vafs.merged = merge(vafs.merged, vafs.merged.cn2, by.x = c(1:length(vafs.merged)), 
            by.y = c(1:length(vafs.merged)), all.x = TRUE)
        vafs.merged = vafs.merged[order(vafs.merged[, 1], vafs.merged[, 
            2]), ]
        print(paste("found", numClusters, "clusters using", clusterMethod, 
            "in full dimensional data"))
    }
    if (!is.null(annotation)) {
        vafs.merged = merge(vafs.merged, annotation, by.x = c(1, 
            2), by.y = c(1, 2), all.x = TRUE, all.y = FALSE)
    }
    return(new("scObject", clust = clust, densities = densityData, 
        dimensions = dimensions, marginalClust = marginalClust, 
        sampleNames = sampleNames, vafs.1d = vafs.1d, vafs.merged = vafs.merged, 
        purities = purities))
}

cleanAndAddCN <- function (vafs, cn, num, cnCallsAreLog2, regionsToExclude, useSexChrs, 
    minimumDepth, copyNumberMargins) {
    names(vafs) = c("chr", "st", "ref", "var", "vaf")
    vafs = vafs[!(vafs$chr == "M" | vafs$chr == "MT"), ]
    if (length(vafs[, 1]) == 0) {
        return(vafs)
    }
    vafs = vafs[!(is.na(vafs$vaf)), ]
    if (length(vafs[, 1]) == 0) {
        return(vafs)
    }
    vafs = unique(vafs)
    for (i in 2:5) {
        if (!(is.numeric(vafs[, i]))) {
            print(paste("ERROR: column", names(vafs)[i], " in sample", 
                i, "is not numeric"))
            stop()
        }
    }
    if (!is.null(regionsToExclude)) {
        vafs = excludeRegions(vafs, regionsToExclude)
        if (length(vafs[, 1]) == 0) {
            return(vafs)
        }
    }
    if (is.null(cn)) {
        vafs$cn = 2
        vafs$cleancn = 2
    }
    else {
        if (cnCallsAreLog2) {
            cn[, 4] = (2^(cn[, 4])) * 2
        }
        vafs = addCnToVafs(vafs, cn, copyNumberMargins)
    }
    vafs = vafs[vafs$vaf > 0 | (vafs$var + vafs$ref) > 0, ]
    vafs$depth = mapply(function(var, ref, vaf) ifelse(vaf == 
        0, var + ref, round(var/(vaf/100))), vafs$var, vafs$ref, 
        vafs$vaf)
    if (!(useSexChrs)) {
        vafs = vafs[vafs$chr != "X" & vafs$chr != "Y", ]
    }
    return(vafs)
}

getDensity <- function (vafs, copyNumberMargins, minimumDepth) {
    densities = vector("list", 4)
    factors = vector("list", 4)
    peakPos = vector("list", 4)
    peakHeights = vector("list", 4)
    maxDensity = 0
    maxDepth = 0
    for (i in 1:4) {
        v = vafs[(vafs$cn > (i - copyNumberMargins)) & (vafs$cn < 
            (i + copyNumberMargins)) & (vafs$depth > minimumDepth), 
            ]
        if (length(v[, 1]) > 0) {
            if (length(v[, 1]) > 1) {
                densities[[i]] = density(v$vaf, from = 0, to = 100, 
                  na.rm = TRUE, adj = 0.85)
                factors[[i]] = (length(v[, 1])/length(vafs[, 
                  1])) * densities[[i]]$y
                p = c(getPeaks(factors[[i]]), FALSE, FALSE)
                peakPos[[i]] = densities[[i]]$x[p]
                peakHeights[[i]] = factors[[i]][p]
                if (max(factors[[i]]) > maxDensity) {
                  maxDensity = max(factors[[i]])
                }
                keep = which(peakHeights[[i]] > maxDensity * 
                  0.1)
                peakPos[[i]] = peakPos[[i]][keep]
                peakHeights[[i]] = peakPos[[i]][keep]
            }
            if (max(v$depth) > maxDepth) {
                maxDepth = max(v$depth)
            }
        }
    }
    return(list(densities = densities, factors = factors, peakPos = peakPos, 
        peakHeights = peakHeights, maxDensity = maxDensity, maxDepth = maxDepth))
}

getPeaks <- function (series, span = 3) {
    z <- embed(series, span)
    s <- span%/%2
    v <- max.col(z) == 1 + s
    result <- c(rep(FALSE, s), v)
    result <- result[1:(length(result) - s)]
    return(result)
}

clusterVafs <- function (vafs.merged, vafMatrix, varMatrix, refMatrix, maximumClusters, 
    method = "bmm", params = NULL, samples = 1, plotIntermediateResults = 0, 
    verbose = 0) {
    apply.uncertainty.self.overlap.condition <- TRUE
    apply.min.items.condition <- TRUE
    apply.overlapping.std.dev.condition <- TRUE
    if (!is.null(params)) {
        if (grepl("no.uncertainty.overlap.detection", params)) {
            cat("Disable uncertainty overlap detection\n")
            apply.uncertainty.self.overlap.condition <- FALSE
        }
        if (grepl("no.min.items.condition", params)) {
            cat("Disable min cluster size detection\n")
            apply.min.items.condition <- FALSE
        }
        if (grepl("no.apply.overlapping.std.dev.condition", params)) {
            cat("Disable overlapping std dev condition\n")
            apply.overlapping.std.dev.condition <- FALSE
        }
    }
    if (method == "bmm") {
        if (!is.null(params)) {
            params = strsplit(params, ", *", perl = TRUE)[[1]]
            overlap.threshold = 0.7
            apply.pvalue.outlier.condition <- TRUE
            outlier.pvalue.threshold = 0.01
            alpha0 = 0.005
            mu0 = 1
            c0 = NULL
            if (grepl("overlap.threshold", params)) {
                overlap.threshold = as.numeric(strsplit(params[grep("overlap.threshold", 
                  params)], " *= *", perl = TRUE)[[1]][2])
            }
            if (grepl("alpha0", params)) {
                alpha0 = as.numeric(strsplit(params[grep("alpha0", 
                  params)], " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using alpha0 = ", alpha0, "\n", sep = ""))
            }
            if (grepl("mu0", params)) {
                mu0 = as.numeric(strsplit(params[grep("mu0", 
                  params)], " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using mu0 = ", mu0, "\n", sep = ""))
            }
            if (grepl("no.pvalue.outlier.detection", params)) {
                cat("Disable pvalue outlier condition\n")
                apply.pvalue.outlier.condition <- FALSE
            }
            if (grepl("outlier.pvalue.threshold", params)) {
                outlier.pvalue.threshold = as.numeric(strsplit(params[grep(" outlier.pvalue.threshold", 
                  params)], " *= *", perl = TRUE)[[1]][2])
            }
            if (grepl("c0", params)) {
                c0 = as.numeric(strsplit(params[grep("c0", params)], 
                  " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using c0 = ", c0, "\n", sep = ""))
            }
            return(clusterWithBmm(vafs.merged, vafMatrix, varMatrix, 
                refMatrix, samples = samples, plotIntermediateResults = plotIntermediateResults, 
                verbose = 0, overlap.threshold = overlap.threshold, 
                apply.pvalue.outlier.condition = apply.pvalue.outlier.condition, 
                initialClusters = maximumClusters, outlier.pvalue.threshold = outlier.pvalue.threshold, 
                apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
                apply.min.items.condition = apply.min.items.condition, 
                apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition, 
                alpha0 = alpha0, mu0 = mu0, c0 = c0))
        }
        return(clusterWithBmm(vafs.merged, vafMatrix, varMatrix, 
            refMatrix, samples = samples, plotIntermediateResults = plotIntermediateResults, 
            verbose = 0, initialClusters = maximumClusters, apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
            apply.min.items.condition = apply.min.items.condition, 
            apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition))
    }
    else if (method == "binomial.bmm") {
        if (!is.null(params)) {
            params = strsplit(params, ", *", perl = TRUE)[[1]]
            if (grepl("overlap.threshold", params)) {
                val = as.numeric(strsplit(params[grep("overlap.threshold", 
                  params)], " *= *", perl = TRUE)[[1]][2])
                return(clusterWithBinomialBmm(vafs.merged, vafMatrix, 
                  varMatrix, refMatrix, samples = samples, plotIntermediateResults = plotIntermediateResults, 
                  verbose = 0, overlap.threshold = val, initialClusters = maximumClusters, 
                  apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
                  apply.min.items.condition = apply.min.items.condition, 
                  apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition))
            }
        }
        return(clusterWithBinomialBmm(vafs.merged, vafMatrix, 
            varMatrix, refMatrix, samples = samples, plotIntermediateResults = plotIntermediateResults, 
            verbose = 0, initialClusters = maximumClusters, apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
            apply.min.items.condition = apply.min.items.condition, 
            apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition))
    }
    else if (method == "gaussian.bmm") {
        W0 = NULL
        beta0 = NULL
        nu0 = NULL
        c0 = NULL
        if (!is.null(params)) {
            if (grepl("W0", params)) {
                W0 = as.numeric(strsplit(params[grep("W0", params)], 
                  " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using W0 = ", W0, "\n", sep = ""))
            }
            if (grepl("beta0", params)) {
                beta0 = as.numeric(strsplit(params[grep("beta0", 
                  params)], " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using beta0 = ", beta0, "\n", sep = ""))
            }
            if (grepl("nu0", params)) {
                nu0 = as.numeric(strsplit(params[grep("nu0", 
                  params)], " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using nu0 = ", nu0, "\n", sep = ""))
            }
            if (grepl("c0", params)) {
                c0 = as.numeric(strsplit(params[grep("c0", params)], 
                  " *= *", perl = TRUE)[[1]][2])
                cat(paste("Using c0 = ", c0, "\n", sep = ""))
            }
        }
        if (!is.null(params)) {
            if (grepl("overlap.threshold", params)) {
                val = as.numeric(strsplit(params[grep("overlap.threshold", 
                  params)], " *= *", perl = TRUE)[[1]][2])
                return(clusterWithGaussianBmm(vafs.merged, vafMatrix, 
                  varMatrix, refMatrix, samples = samples, plotIntermediateResults = plotIntermediateResults, 
                  verbose = 0, overlap.threshold = val, initialClusters = maximumClusters, 
                  apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
                  apply.min.items.condition = apply.min.items.condition, 
                  apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition, 
                  nu0 = nu0, beta0 = beta0, W0 = W0, c0 = c0))
            }
        }
        return(clusterWithGaussianBmm(vafs.merged, vafMatrix, 
            varMatrix, refMatrix, samples = samples, plotIntermediateResults = plotIntermediateResults, 
            verbose = 0, initialClusters = maximumClusters, apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
            apply.min.items.condition = apply.min.items.condition, 
            apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition, 
            nu0 = nu0, beta0 = beta0, W0 = W0, c0 = c0))
    }
    else {
        print("Error: please choose a supported clustering method\n[bmm|gaussian.bmm|binomial.bmm]")
        return(0)
    }
}

clusterWithBmm <- function (vafs.merged, vafs, vars, refs, initialClusters = 10, 
    samples = 1, plotIntermediateResults = 0, verbose = TRUE, 
    overlap.threshold = 0.7, apply.pvalue.outlier.condition = TRUE, 
    outlier.pvalue.threshold = 0.01, apply.uncertainty.self.overlap.condition = TRUE, 
    apply.min.items.condition = TRUE, apply.overlapping.std.dev.condition = TRUE, 
    alpha0 = 0.005, mu0 = 1, c0 = NULL) {
    suppressPackageStartupMessages(library(bmm))
    initialClusters = initialClusters
    if (length(vafs[, 1]) <= initialClusters) {
        print(paste("ERROR: only", length(vafs[, 1]), " points - not enough points to cluster when using", 
            initialClusters, "intialClusters. Provide more data or reduce your maximumClusters option"))
        return(list(NULL))
    }
    shift.by.machine.precision <- TRUE
    delta <- 100 * .Machine$double.eps
    if (shift.by.machine.precision) {
        vafs[which(vafs <= delta)] = delta
    }
    else {
        num.dimensions <- ncol(vafs)
        for (m in 1:num.dimensions) {
            flag <- vafs[, m] <= delta
            min.vaf <- min(vafs[!flag, m])
            vafs[flag, m] <- unlist(lapply(refs[flag, m], function(x) min(min.vaf, 
                1/x)))
        }
    }
    vafs[which(vafs >= (1 - delta))] = 1 - delta
    hyperparams <- init.bmm.hyperparameters(vafs, initialClusters)
    if (!is.null(c0)) {
        hyperparams$c0 <- rep(c0, length(hyperparams$c0))
    }
    hyperparams$alpha0 <- matrix(data = alpha0, nrow = nrow(hyperparams$alpha0), 
        ncol = ncol(hyperparams$alpha0))
    hyperparams$beta0 <- hyperparams$alpha0
    hyperparams$mu0 <- matrix(data = mu0, nrow = nrow(hyperparams$mu0), 
        ncol = ncol(hyperparams$mu0))
    hyperparams$nu0 <- hyperparams$mu0
    params <- init.bmm.parameters(vafs, initialClusters, hyperparams$mu0, 
        hyperparams$alpha0, hyperparams$nu0, hyperparams$beta0, 
        hyperparams$c0, verbose = TRUE)
    bmm.results <- bmm.filter.clusters(vafs.merged, vafs, initialClusters, 
        params$r, params$mu, params$alpha, params$nu, params$beta, 
        params$c, params$E.pi, hyperparams$mu0, hyperparams$alpha0, 
        hyperparams$nu0, hyperparams$beta0, hyperparams$c0, convergence.threshold = 10^-4, 
        max.iterations = 10000, verbose = verbose, plotIntermediateResults = plotIntermediateResults, 
        overlap.threshold = overlap.threshold, apply.pvalue.outlier.condition = apply.pvalue.outlier.condition, 
        outlier.pvalue.threshold = outlier.pvalue.threshold, 
        apply.uncertainty.self.overlap.condition = apply.uncertainty.self.overlap.condition, 
        apply.min.items.condition = apply.min.items.condition, 
        apply.overlapping.std.dev.condition = apply.overlapping.std.dev.condition)
    if (bmm.results$retVal != 0) {
        cat("WARNING: bmm failed to converge. No clusters assigned\n")
        return(NULL)
    }
    probs = bmm.results$r
    numPoints = length(probs[, 1])
    numClusters = length(probs[1, ])
    clusters = hardClusterAssignments(numPoints, 1:numClusters, 
        probs)
    intervals = bmm.narrowest.mean.interval.about.centers(bmm.results$mu, 
        bmm.results$alpha, bmm.results$nu, bmm.results$beta, 
        0.68)
    means = intervals$centers
    if (verbose) {
        print("Cluster Centers:")
        print(means)
    }
    lower = intervals$lb
    upper = intervals$ub
    if (dim(bmm.results$outliers)[1] > 0) {
        cat("Outliers:\n")
        print(bmm.results$outliers)
    }
    n <- 1000
    x <- seq(0, 1, 1/n)[2:n]
    n = n - 1
    y <- rep.int(0, n)
    y = t(matrix(rep(y, dim(vafs)[2]), ncol = dim(vafs)[2]))
    yms = list()
    for (k in 1:numClusters) {
        yms[[k]] <- y
    }
    for (dim in 1:dim(vafs)[2]) {
        ym <- matrix(data = 0, nrow = numClusters, ncol = n)
        num.iterations <- 100
        for (k in 1:numClusters) {
            for (i in 1:n) {
                ym[k, i] <- bmm.component.posterior.predictive.density(x[i], 
                  bmm.results$mu[dim, k], bmm.results$alpha[dim, 
                    k], bmm.results$nu[dim, k], bmm.results$beta[dim, 
                    k], bmm.results$E.pi[k], num.samples = num.iterations)
                yms[[k]][dim, i] <- ym[k, i]
                y[dim, i] <- y[dim, i] + ym[k, i]
            }
        }
    }
    x = x * 100
    cat("Converged on the following parameters:\n")
    cat("mu:\n")
    write.table(bmm.results$mu, col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    cat("alpha:\n")
    write.table(bmm.results$alpha, col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    cat("nu:\n")
    write.table(bmm.results$nu, col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    cat("beta:\n")
    write.table(bmm.results$beta, col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    cat("pi:\n")
    cat(paste(bmm.results$E.pi, collapse = "\t"))
    cat("\n")
    return(list(cluster.assignments = clusters, cluster.probabilities = probs, 
        cluster.means = means, cluster.upper = upper, cluster.lower = lower, 
        fit.x = x, fit.y = y, individual.fits.y = yms, mu = bmm.results$mu, 
        alpha = bmm.results$alpha, nu = bmm.results$nu, beta = bmm.results$beta, 
        pi = bmm.results$E.pi, cluster.method = "bmm"))
}

bmm.filter.clusters <- function (vafs.merged, X, N.c, r, mu, alpha, nu, beta, c, E.pi, 
    mu0, alpha0, nu0, beta0, c0, convergence.threshold = 10^-4, 
    max.iterations = 10000, verbose = 0, plotIntermediateResults = 0, 
    overlap.threshold = 0.7, apply.pvalue.outlier.condition = TRUE, 
    outlier.pvalue.threshold = 0.01, apply.uncertainty.self.overlap.condition = TRUE, 
    apply.min.items.condition = TRUE, apply.overlapping.std.dev.condition = TRUE) {
    total.iterations <- 0
    N <- dim(X)[1]
    num.dimensions <- dim(X)[2]
    effective.overlap.threshold = (overlap.threshold^(1/num.dimensions))
    cat("Using threshold: ", effective.overlap.threshold, "\n")
    cluster.names <- 1:N.c
    x.colnames <- colnames(X)
    outliers <- matrix(data = 0, nrow = 0, ncol = dim(X)[2])
    colnames(outliers) <- colnames(X)
    E.pi.prev <- rep(0, N.c)
    E.pi.prev <- E.pi
    width <- as.double(erf(1.5/sqrt(2)))
    if (plotIntermediateResults > 0) {
        probs <- r
        numClusters = length(probs[1, ])
        numPoints = length(probs[, 1])
        clusters = hardClusterAssignments(numPoints, cluster.names, 
            probs)
        ellipse.width <- as.double(erf(1/sqrt(2)))
        SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, 
            alpha, nu, beta, ellipse.width)
        SEM.centers <- 100 * t(SEM.res$centers)
        SEMs.lb <- 100 * t(SEM.res$lb)
        SEMs.ub <- 100 * t(SEM.res$ub)
        std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
            alpha, nu, beta, ellipse.width)
        std.dev.centers <- 100 * t(std.dev.res$centers)
        std.dev.lb <- 100 * t(std.dev.res$lb)
        std.dev.ub <- 100 * t(std.dev.res$ub)
        ellipse.metadata <- list(SEMs.lb = SEMs.lb, SEMs.ub = SEMs.ub, 
            std.dev.lb = std.dev.lb, std.dev.ub = std.dev.ub)
        means = SEM.res$centers
        clust <- list(cluster.assignments = clusters, cluster.probabilities = probs, 
            cluster.means = means, cluster.upper = means, cluster.lower = means, 
            fit.x = NULL, fit.y = NULL, individual.fits.y = NULL)
        vafs.with.assignments = cbind(vafs.merged, cluster = clust$cluster.assignments)
        vafs.with.assignments = cbind(vafs.with.assignments, 
            cluster.prob = clust$cluster.probabilities)
        outputFile <- paste("iter.", total.iterations, ".pdf", 
            sep = "")
        sampleNames <- names(vafs.merged)
        sampleNames <- sampleNames[grepl(pattern = ".ref", sampleNames)]
        for (s in 1:length(sampleNames)) {
            sampleNames[s] <- substr(sampleNames[s], 1, nchar(sampleNames[s]) - 
                4)
        }
        positionsToHighlight <- NULL
        highlightsHaveNames <- FALSE
        overlayClusters <- TRUE
        sco <- new("scObject", dimensions = num.dimensions, sampleNames = sampleNames, 
            vafs.merged = vafs.with.assignments)
        sc.plot2d(sco, outputFile, ellipse.metadata = ellipse.metadata, 
            positionsToHighlight = positionsToHighlight, highlightsHaveNames = highlightsHaveNames)
    }
    while (TRUE) {
        if (plotIntermediateResults > 0) {
            max.iterations <- plotIntermediateResults
        }
        bmm.res <- bmm.fixed.num.components(X, N.c, r, mu, alpha, 
            nu, beta, c, E.pi, mu0, alpha0, nu0, beta0, c0, convergence.threshold, 
            max.iterations, verbose)
        if ((bmm.res$retVal != 0) & (plotIntermediateResults == 
            0)) {
            stop("Failed to converge!\n")
        }
        mu <- bmm.res$mu
        alpha <- bmm.res$alpha
        nu <- bmm.res$nu
        beta <- bmm.res$beta
        c <- bmm.res$c
        E.pi <- bmm.res$E.pi
        ubar <- bmm.res$ubar
        vbar <- bmm.res$vbar
        r <- bmm.res$r
        total.iterations <- total.iterations + bmm.res$num.iterations
        ln.rho <- bmm.res$ln.rho
        E.lnu <- bmm.res$E.lnu
        E.lnv <- bmm.res$E.lnv
        E.lnpi <- bmm.res$E.lnpi
        E.quadratic.u <- bmm.res$E.quadratic.u
        E.quadratic.v <- bmm.res$E.quadratic.v
        if (plotIntermediateResults > 0) {
            probs <- r
            numPoints = length(probs[, 1])
            numClusters = length(probs[1, ])
            clusters = hardClusterAssignments(numPoints, cluster.names, 
                probs)
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, 
                alpha, nu, beta, ellipse.width)
            SEM.centers <- 100 * t(SEM.res$centers)
            SEMs.lb <- 100 * t(SEM.res$lb)
            SEMs.ub <- 100 * t(SEM.res$ub)
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
                alpha, nu, beta, ellipse.width)
            std.dev.centers <- 100 * t(std.dev.res$centers)
            std.dev.lb <- 100 * t(std.dev.res$lb)
            std.dev.ub <- 100 * t(std.dev.res$ub)
            ellipse.metadata <- list(SEMs.lb = SEMs.lb, SEMs.ub = SEMs.ub, 
                std.dev.lb = std.dev.lb, std.dev.ub = std.dev.ub)
            means = SEM.res$centers
            clust <- list(cluster.assignments = clusters, cluster.probabilities = probs, 
                cluster.means = means, cluster.upper = means, 
                cluster.lower = means, fit.x = NULL, fit.y = NULL, 
                individual.fits.y = NULL)
            vafs.with.assignments = cbind(vafs.merged, cluster = clust$cluster.assignments)
            vafs.with.assignments = cbind(vafs.with.assignments, 
                cluster.prob = clust$cluster.probabilities)
            outputFile <- paste("iter.", total.iterations, ".pdf", 
                sep = "")
            positionsToHighlight <- NULL
            highlightsHaveNames <- FALSE
            overlayClusters <- TRUE
            sco <- new("scObject", dimensions = num.dimensions, 
                sampleNames = sampleNames, vafs.merged = vafs.with.assignments)
            sc.plot2d(sco, outputFile, ellipse.metadata = ellipse.metadata, 
                positionsToHighlight = positionsToHighlight, 
                highlightsHaveNames = highlightsHaveNames)
        }
        do.inner.iteration <- FALSE
        apply.large.SEM.condition <- FALSE
        apply.overlapping.SEM.condition <- FALSE
        apply.outlier.condition <- FALSE
        indices.to.keep <- rep(TRUE, N.c)
        remove.data <- FALSE
        remove.data.from.small.clusters <- FALSE
        clusters <- hardClusterAssignments(N, 1:N.c, r)
        if ((apply.min.items.condition == TRUE) & (N.c > 1)) {
            threshold.pts <- max(3, ceiling(0.005 * N))
            threshold.pts <- 3
            num.items.per.cluster <- rep(0, N.c)
            for (n in 1:N) {
                num.items.per.cluster[clusters[n]] <- num.items.per.cluster[clusters[n]] + 
                  1
            }
            indices.above.threshold <- num.items.per.cluster >= 
                threshold.pts
            if (any(indices.above.threshold == FALSE)) {
                dropped.clusters <- (1:N.c)[!indices.above.threshold]
                expected.means <- ubar/(ubar + vbar)
                save.verbose <- verbose
                verbose <- TRUE
                for (k in dropped.clusters) {
                  indices.to.keep[k] <- FALSE
                  if (verbose) {
                    cat("Dropped cluster", k, "with too few variants (", 
                      num.items.per.cluster[k], ") center: ")
                    for (n in 1:length(expected.means[, k])) {
                      cat(expected.means[n, k])
                      if (length(expected.means[, k]) > 1) {
                        cat(", ")
                      }
                    }
                    cat("\n")
                  }
                  break
                }
                verbose <- save.verbose
            }
            if (remove.data.from.small.clusters == TRUE) {
                if (any(indices.to.keep == FALSE)) {
                  overlaps <- calculate.self.overlap(r)
                  if (any(!is.na(overlaps[!indices.to.keep]) & 
                    (overlaps[!indices.to.keep] > 0.99))) {
                    remove.data <- TRUE
                    indices.to.keep <- indices.to.keep | (is.na(overlaps) | 
                      (overlaps < 0.99))
                  }
                }
            }
        }
        if (all(indices.to.keep == TRUE) & (apply.outlier.condition == 
            TRUE)) {
            ellipse.width <- as.double(erf(1.5/sqrt(2)))
            ellipse.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
                alpha, nu, beta, ellipse.width)
            ellipse.centers <- t(ellipse.res$centers)
            ellipse.lb <- t(ellipse.res$lb)
            ellipse.ub <- t(ellipse.res$ub)
            removed.pt <- FALSE
            for (k in 1:N.c) {
                if (verbose) {
                  cat("Cluster", k, "\n")
                  print(ellipse.lb[k, ])
                  print(ellipse.ub[k, ])
                }
                indices <- (1:N)[clusters == k]
                remove.i <- TRUE
                for (i in indices) {
                  for (m in 1:num.dimensions) {
                    lower <- ellipse.lb[k, m]
                    upper <- ellipse.ub[k, m]
                    if (!(X[i, m] < lower) & !(X[i, m] > upper)) {
                      remove.i <- FALSE
                      break
                    }
                    if ((X[i, m] < lower) | (X[i, m] > upper)) {
                      if (verbose) {
                        cat("Point ")
                        for (n in 1:num.dimensions) {
                          cat(X[i, n], " ")
                        }
                        cat("is outside of range.\n")
                        print(ellipse.lb[k, ])
                        print(ellipse.ub[k, ])
                      }
                      removed.pt <- TRUE
                      do.inner.iteration <- TRUE
                      outliers <- rbind(outliers, X[i, , drop = F])
                      X[i, ] <- NA
                      break
                    }
                  }
                  if (remove.i) {
                    removed.pt <- TRUE
                    do.inner.iteration <- TRUE
                    if (verbose) {
                      cat("Point ")
                      for (n in 1:num.dimensions) {
                        cat(X[i, n], " ")
                      }
                      cat("is outside of range.\n")
                      print(ellipse.lb[k, ])
                      print(ellipse.ub[k, ])
                    }
                    outliers <- rbind(outliers, X[i, , drop = F])
                    X[i, ] <- NA
                  }
                }
            }
            if (removed.pt == TRUE) {
                next
            }
        }
        if (all(indices.to.keep == TRUE) & (apply.pvalue.outlier.condition == 
            TRUE)) {
            ellipse.width <- as.double(erf(0.75/sqrt(2)))
            ellipse.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
                alpha, nu, beta, ellipse.width)
            ellipse.centers <- t(ellipse.res$centers)
            ellipse.lb <- t(ellipse.res$lb)
            ellipse.ub <- t(ellipse.res$ub)
            removed.pt <- FALSE
            for (k in 1:N.c) {
                scrutinize.i <- TRUE
                indices <- (1:N)[clusters == k]
                for (i in indices) {
                  for (m in 1:num.dimensions) {
                    lower <- ellipse.lb[k, m]
                    upper <- ellipse.ub[k, m]
                    if (!(X[i, m] < lower) & !(X[i, m] > upper)) {
                      scrutinize.i <- FALSE
                      break
                    }
                  }
                  if (scrutinize.i) {
                    num.samples <- 1000
                    proportions <- matrix(data = 0, nrow = num.samples, 
                      ncol = num.dimensions)
                    for (m in 1:num.dimensions) {
                      proportions[, m] <- sample.bmm.component.proportion(mu[m, 
                        k], alpha[m, k], nu[m, k], beta[m, k], 
                        num.samples)
                    }
                    probabilities <- rep(1, num.samples)
                    for (n in 1:num.samples) {
                      for (m in 1:num.dimensions) {
                        probabilities[n] <- probabilities[n] * 
                          bmm.component.posterior.predictive.density(proportions[n, 
                            m], mu[m, k], alpha[m, k], nu[m, 
                            k], beta[m, k], num.samples)
                      }
                    }
                    i.prob <- 1
                    for (m in 1:num.dimensions) {
                      i.prob <- i.prob * bmm.component.posterior.predictive.density(X[i, 
                        m], mu[m, k], alpha[m, k], nu[m, k], 
                        beta[m, k], num.samples)
                    }
                    pvalue <- length((1:length(probabilities))[probabilities <= 
                      i.prob])/length(probabilities)
                    if (verbose) {
                      cat("Point (pvalue = ", pvalue, ") ", sep = "")
                      for (n in 1:num.dimensions) {
                        cat(X[i, n], " ")
                      }
                      cat("\n")
                    }
                    if (pvalue < outlier.pvalue.threshold) {
                      removed.pt <- TRUE
                      do.inner.iteration <- TRUE
                      if (verbose) {
                        cat("Point ")
                        for (n in 1:num.dimensions) {
                          cat(X[i, n], " ")
                        }
                        cat(" has small pvalue = ", pvalue, "\n")
                      }
                      outliers <- rbind(outliers, X[i, , drop = F])
                      X[i, ] <- NA
                    }
                  }
                }
            }
            if (removed.pt == TRUE) {
                next
            }
        }
        if (all(indices.to.keep == TRUE) & (apply.uncertainty.self.overlap.condition == 
            TRUE) & (N.c > 1)) {
            overlaps <- calculate.self.overlap(r)
            if (min(overlaps, na.rm = TRUE) < effective.overlap.threshold) {
                for (k in 1:N.c) {
                  if (is.nan(overlaps[k])) {
                    next
                  }
                  if ((overlaps[k] < effective.overlap.threshold) & 
                    (overlaps[k] == min(overlaps, na.rm = TRUE))) {
                    indices.to.keep[k] <- FALSE
                    if (verbose) {
                      cat(sprintf("Condition (%dD): Remove cluster %d pi = %.3f self-overlap = %.3f", 
                        num.dimensions, k, E.pi[k], overlaps[k]))
                    }
                    expected.means <- ubar/(ubar + vbar)
                    if (verbose) {
                      cat(" center: ")
                      for (n in 1:length(expected.means[, k])) {
                        cat(expected.means[n, k])
                        if (length(expected.means[, k]) > 1) {
                          cat(", ")
                        }
                      }
                      cat("\n")
                    }
                    break
                  }
                }
            }
            if (verbose) {
                for (k in 1:N.c) {
                  cat(sprintf("Cluster %d pi = %.3f self-overlap = %.3f\n", 
                    k, E.pi[k], overlaps[k]))
                }
            }
        }
        if (all(indices.to.keep == TRUE) & (apply.large.SEM.condition == 
            TRUE) & (N.c > 1)) {
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, 
                alpha, nu, beta, width)
            SEM.centers <- t(SEM.res$centers)
            SEMs.lb <- t(SEM.res$lb)
            SEMs.ub <- t(SEM.res$ub)
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
                alpha, nu, beta, width)
            std.dev.centers <- t(std.dev.res$centers)
            std.dev.lb <- t(std.dev.res$lb)
            std.dev.ub <- t(std.dev.res$ub)
            for (k in 1:N.c) {
                if (verbose) {
                  cat(sprintf("%dD: Cluster %d pi = %.3f: ", 
                    num.dimensions, k, E.pi[k]))
                }
                greater.than.30 <- TRUE
                greater.than.02 <- TRUE
                SEM.width.threshold <- 0.02
                for (m in 1:num.dimensions) {
                  center <- SEM.centers[k, m]
                  lower <- SEMs.lb[k, m]
                  upper <- SEMs.ub[k, m]
                  std.dev.width <- (std.dev.ub[k, m] - std.dev.lb[k, 
                    m])
                  SEM.width <- (SEMs.ub[k, m] - SEMs.lb[k, m])
                  if ((SEM.width/std.dev.width) <= 0.3) {
                    greater.than.30 <- FALSE
                  }
                  if (SEM.width <= SEM.width.threshold) {
                    greater.than.02 <- FALSE
                  }
                  if (verbose) {
                    cat(sprintf("%.3f (%.3f, %.3f) [(%.3f) %.3f, %.3f] {%.3f} ", 
                      center, lower, upper, width, std.dev.lb[k, 
                        m], std.dev.ub[k, m], SEM.width/std.dev.width))
                    if (greater.than.30) {
                      cat("* ")
                    }
                    if (greater.than.02) {
                      cat("**")
                    }
                  }
                }
                if (verbose) {
                  cat("\n")
                }
                if (greater.than.30 & greater.than.02) {
                  indices.to.keep[k] <- FALSE
                }
            }
            if (any(indices.to.keep == FALSE)) {
                remove.data <- TRUE
            }
        }
        if (all(indices.to.keep == TRUE) & (apply.overlapping.SEM.condition == 
            TRUE) & (N.c > 1)) {
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, 
                alpha, nu, beta, width)
            SEM.centers <- t(SEM.res$centers)
            SEMs.lb <- t(SEM.res$lb)
            SEMs.ub <- t(SEM.res$ub)
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
                alpha, nu, beta, width)
            std.dev.centers <- t(std.dev.res$centers)
            std.dev.lb <- t(std.dev.res$lb)
            std.dev.ub <- t(std.dev.res$ub)
            pi.threshold = 10^-2
            for (i in 1:N.c) {
                i.subsumed.by.another.cluster <- FALSE
                for (i2 in 1:N.c) {
                  if (i == i2) {
                    next
                  }
                  if (indices.to.keep[i2] == FALSE) {
                    next
                  }
                  if (E.pi[i2] < pi.threshold) {
                    next
                  }
                  i.subsumed.by.another.cluster <- TRUE
                  for (l in 1:num.dimensions) {
                    i.center <- std.dev.centers[i, l]
                    if ((i.center < std.dev.lb[i2, l]) | (i.center > 
                      std.dev.ub[i2, l])) {
                      i.subsumed.by.another.cluster <- FALSE
                      break
                    }
                  }
                  if (i.subsumed.by.another.cluster == TRUE) {
                    if (verbose) {
                      cat(sprintf("2. Dropping cluster with center: "))
                      for (l in 1:num.dimensions) {
                        cat(sprintf("%.3f ", std.dev.centers[i, 
                          l]))
                      }
                      cat(sprintf("because it overlaps with: "))
                      for (l in 1:num.dimensions) {
                        cat(sprintf("(%.3f, %.3f) ", std.dev.lb[i2, 
                          l], std.dev.ub[i2, l]))
                      }
                      cat("\n")
                      break
                    }
                  }
                }
                if (i.subsumed.by.another.cluster == TRUE) {
                  indices.to.keep[i] <- FALSE
                }
            }
        }
        if (all(indices.to.keep == TRUE) & (apply.overlapping.std.dev.condition == 
            TRUE) & (N.c > 1)) {
            SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, 
                alpha, nu, beta, width)
            SEM.centers <- t(SEM.res$centers)
            SEMs.lb <- t(SEM.res$lb)
            SEMs.ub <- t(SEM.res$ub)
            std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
                alpha, nu, beta, width)
            std.dev.centers <- t(std.dev.res$centers)
            std.dev.lb <- t(std.dev.res$lb)
            std.dev.ub <- t(std.dev.res$ub)
            overlaps <- matrix(data = 1, nrow = N.c, ncol = N.c)
            std.dev.overlap.threshold <- 0
            for (i in 1:N.c) {
                for (i2 in 1:N.c) {
                  if (verbose) {
                    cat("   ")
                  }
                  for (l in 1:num.dimensions) {
                    fraction <- 0
                    if ((std.dev.lb[i, l] < std.dev.ub[i2, l]) & 
                      (std.dev.ub[i, l] > std.dev.lb[i2, l])) {
                      overlap <- min(std.dev.ub[i, l], std.dev.ub[i2, 
                        l]) - max(std.dev.lb[i, l], std.dev.lb[i2, 
                        l])
                      fraction <- overlap/(std.dev.ub[i, l] - 
                        std.dev.lb[i, l])
                    }
                    if ((fraction > std.dev.overlap.threshold) & 
                      (overlaps[i, i2] != 0)) {
                      overlaps[i, i2] <- fraction
                    }
                    else {
                      overlaps[i, i2] <- 0
                    }
                  }
                }
            }
            save.verbose <- verbose
            verbose <- TRUE
            for (i in 2:N.c) {
                if (indices.to.keep[i] == FALSE) {
                  next
                }
                for (i2 in 1:(i - 1)) {
                  if (indices.to.keep[i2] == FALSE) {
                    next
                  }
                  if (overlaps[i, i2] > 0) {
                    if ((overlaps[i2, i] > 0) & (E.pi[i2] < E.pi[i])) {
                      if (verbose) {
                        cat("Condition (", num.dimensions, "D): Removing ", 
                          i2, " because of overlap (", overlaps[i2, 
                            i], ") with i = ", i, "\n")
                      }
                      indices.to.keep[i2] <- FALSE
                    }
                    else {
                      if (verbose) {
                        cat("Condition (", num.dimensions, "D): Removing ", 
                          i, " because of overlap (", overlaps[i2, 
                            i], ") with i2 = ", i2, "\n")
                      }
                      indices.to.keep[i] <- FALSE
                    }
                  }
                }
            }
            verbose <- save.verbose
        }
        if (any(indices.to.keep == FALSE)) {
            do.inner.iteration <- TRUE
            numeric.indices <- (1:N.c)
            cluster.names <- cluster.names[indices.to.keep]
            E.pi <- E.pi[indices.to.keep]
            E.lnpi <- E.lnpi[indices.to.keep]
            if (remove.data == TRUE) {
                cluster.indices.to.keep <- (1:N.c)[indices.to.keep]
                new.outliers <- matrix(X[!(clusters == 0) & !(clusters %in% 
                  cluster.indices.to.keep), ], ncol = num.dimensions)
                outliers <- rbind(outliers, new.outliers)
                X[!(clusters == 0) & !(clusters %in% cluster.indices.to.keep), 
                  ] <- NA
            }
            colnames(X) <- x.colnames
            N.c <- length(E.pi)
            E.pi.prev <- E.pi.prev[indices.to.keep]
            c <- c[indices.to.keep]
            c0 <- c0[indices.to.keep]
            r <- matrix(r[, indices.to.keep], nrow = N, ncol = N.c)
            ln.rho <- matrix(ln.rho[, indices.to.keep], nrow = N, 
                ncol = N.c)
            for (n in 1:N) {
                if (any(is.na(ln.rho[n, ]))) {
                  r[n, ] <- rep(NA, N.c)
                  next
                }
                row.sum <- log(sum(exp(ln.rho[n, ] - max(ln.rho[n, 
                  ])))) + max(ln.rho[n, ])
                for (k in 1:N.c) {
                  r[n, k] = exp(ln.rho[n, k] - row.sum)
                }
            }
            mu <- matrix(mu[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            nu <- matrix(nu[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            mu0 <- matrix(mu0[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            nu0 <- matrix(nu0[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            alpha <- matrix(alpha[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            beta <- matrix(beta[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            alpha0 <- matrix(alpha0[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            beta0 <- matrix(beta0[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            ubar <- matrix(ubar[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            vbar <- matrix(vbar[, indices.to.keep], nrow = num.dimensions, 
                ncol = N.c)
            E.pi.prev <- E.pi
            E.lnpi <- E.lnpi[indices.to.keep]
            E.lnu <- E.lnu[indices.to.keep]
            E.lnv <- E.lnv[indices.to.keep]
            E.quadratic.u <- matrix(E.quadratic.u[, indices.to.keep], 
                nrow = num.dimensions, ncol = N.c)
            E.quadratic.v <- matrix(E.quadratic.v[, indices.to.keep], 
                nrow = num.dimensions, ncol = N.c)
        }
        if (do.inner.iteration == FALSE) {
            break
        }
    }
    SEM.res <- bmm.narrowest.mean.interval.about.centers(mu, 
        alpha, nu, beta, width)
    SEM.centers <- t(SEM.res$centers)
    SEMs.lb <- t(SEM.res$lb)
    SEMs.ub <- t(SEM.res$ub)
    std.dev.res <- bmm.narrowest.proportion.interval.about.centers(mu, 
        alpha, nu, beta, width)
    std.dev.centers <- t(std.dev.res$centers)
    std.dev.lb <- t(std.dev.res$lb)
    std.dev.ub <- t(std.dev.res$ub)
    save.verbose <- verbose
    verbose <- TRUE
    for (k in 1:N.c) {
        if (verbose) {
            cat(sprintf("Cluster %d pi = %.3f center =", k, E.pi[k]))
        }
        for (d in 1:num.dimensions) {
            center <- SEM.centers[k, d]
            if (verbose) {
                cat(sprintf(" %.3f", center))
            }
        }
        if (verbose) {
            cat(" SEM =")
        }
        for (d in 1:num.dimensions) {
            lb <- SEMs.lb[k, d]
            ub <- SEMs.ub[k, d]
            if (verbose) {
                cat(sprintf(" (%.3f, %.3f)", lb, ub))
            }
        }
        if (verbose) {
            cat(" sd =")
        }
        for (d in 1:num.dimensions) {
            lb <- std.dev.lb[k, d]
            ub <- std.dev.ub[k, d]
            if (verbose) {
                cat(sprintf(" (%.3f, %.3f)", lb, ub))
            }
        }
        if (verbose) {
            cat("\n")
        }
    }
    verbose <- save.verbose
    if (verbose) {
        cat(sprintf("total iterations = %d\n", total.iterations))
    }
    retList <- list(retVal = 0, mu = mu, alpha = alpha, nu = nu, 
        beta = beta, c = c, r = r, num.iterations = total.iterations, 
        ln.rho = ln.rho, E.lnu = E.lnu, E.lnv = E.lnv, E.lnpi = E.lnpi, 
        E.pi = E.pi, E.quadratic.u = E.quadratic.u, E.quadratic.v = E.quadratic.v, 
        ubar = ubar, vbar = vbar, outliers = outliers)
    return(retList)
}

hardClusterAssignments <- function (numPoints, cluster.names, probabilities) {
    assignments <- rep(0, numPoints)
    for (n in 1:numPoints) {
        max.cluster <- 0
        max.assignment <- -1
        for (k in 1:length(cluster.names)) {
            if (!is.na(probabilities[n, k]) & (probabilities[n, 
                k] > max.assignment)) {
                max.assignment <- probabilities[n, k]
                max.cluster <- cluster.names[k]
            }
        }
        assignments[n] <- max.cluster
    }
    return(assignments)
}

calculate.self.overlap <- function (r) {
    N.c <- dim(r)[2]
    N <- dim(r)[1]
    overlaps <- rep(0, N.c)
    ones <- rep(1, N)
    for (k in 1:N.c) {
        overlaps[k] <- sum(r[, k] * r[, k], na.rm = TRUE)/sum(ones * 
            r[, k], na.rm = TRUE)
    }
    return(overlaps)
}

reorderClust <- function (clust) {
    nums = unique(clust$cluster.assignments)[unique(clust$cluster.assignments) > 
        0]
    dist = c()
    for (i in nums) {
        z = clust$cluster.means
        dist = c(dist, sqrt(sum(clust$cluster.means[, i]^2)))
    }
    df = data.frame(nums, dist = dist)
    df = cbind(df[rev(order(df$dist)), ], new = 1:length(nums))
    ass = list()
    prob = list()
    means = list()
    upper = list()
    lower = list()
    individual = list()
    for (i in 1:length(df[, 1])) {
        oldnum = df[i, 1]
        ass[[i]] = which(clust$cluster.assignments == oldnum)
        prob[[i]] = clust$cluster.probabilities[, oldnum]
        means[[i]] = clust$cluster.means[, oldnum]
        lower[[i]] = clust$cluster.lower[, oldnum]
        upper[[i]] = clust$cluster.upper[, oldnum]
        individual[[i]] = clust$individual.fits.y[[oldnum]]
        if (is.null(clust$individual.fits.y[[oldnum]])) {
            individual[[i]] = 0
        }
    }
    for (i in 1:length(df[, 1])) {
        clust$cluster.assignments[ass[[i]]] = i
        clust$cluster.probabilities[, i] = prob[[i]]
        clust$cluster.means[, i] = means[[i]]
        clust$cluster.lower[, i] = lower[[i]]
        clust$cluster.upper[, i] = upper[[i]]
        clust$individual.fits.y[[i]] = individual[[i]]
    }
    if (clust$cluster.method == "bmm") {
        clust <- reorderBetaClust(clust, df)
    }
    else if (clust$cluster.method == "gaussian.bmm") {
        clust <- reorderGaussianClust(clust, df)
    }
    else if (clust$cluster.method == "binomial.bmm") {
        clust <- reorderBinomialClust(clust, df)
    }
    else {
        stop(paste("Cluster reordering not implemented for method", 
            clust$cluster.method))
    }
    return(clust)
}

reorderBetaClust <- function (clust, ord) {
    mu = list()
    alpha = list()
    nu = list()
    beta = list()
    pi = list()
    for (i in 1:length(ord[, 1])) {
        oldnum = ord[i, 1]
        mu[[i]] = clust$mu[, oldnum]
        alpha[[i]] = clust$alpha[, oldnum]
        nu[[i]] = clust$nu[, oldnum]
        beta[[i]] = clust$beta[, oldnum]
        pi[[i]] = clust$pi[oldnum]
    }
    for (i in 1:length(ord[, 1])) {
        clust$mu[, i] = mu[[i]]
        clust$alpha[, i] = alpha[[i]]
        clust$nu[, i] = nu[[i]]
        clust$beta[, i] = beta[[i]]
        clust$pi[i] = pi[[i]]
    }
    return(clust)
}

sc.plot1d <- function (sco, outputFile = NULL, cnToPlot = c(1, 2, 3, 4), showCopyNumberScatterPlots = TRUE, 
    highlightSexChrs = TRUE, positionsToHighlight = NULL, highlightsHaveNames = FALSE, 
    overlayClusters = TRUE, overlayIndividualModels = TRUE, showHistogram = FALSE, 
    showTitle = TRUE, biggerText = FALSE, highlightsOnHistogram = FALSE, 
    highlightCnPoints = FALSE) {
    densityData = sco@densities
    vafs.merged = sco@vafs.merged
    sampleNames = sco@sampleNames
    dimensions = sco@dimensions
    if (max(cnToPlot) > 4 | min(cnToPlot) < 1) {
        print("sciClone supports plotting of copy numbers between 1 and 4 at this time")
    }
    minimumLabelledPeakHeight = 999
    onlyLabelHighestPeak = TRUE
    add.legend <- FALSE
    addpts <- NULL
    if (!is.null(positionsToHighlight)) {
        names(positionsToHighlight) = c("chr", "st", "name")
        addpts = merge(vafs.merged, positionsToHighlight, by.x = c("chr", 
            "st"), by.y = c("chr", "st"))
        if (showCopyNumberScatterPlots & (length(cnToPlot) < 
            2) & highlightsHaveNames) {
            add.legend <- TRUE
        }
    }
    if (highlightsHaveNames & add.legend) {
        if (!(is.null(positionsToHighlight))) {
            if (length(positionsToHighlight) < 3) {
                print("ERROR: named plot requires names in the third column of the positionsToHighlight data frame")
                return(0)
            }
            cnToPlot = c(2)
        }
        else {
            print("ERROR: highlightsHaveNames requires a 3-column dataframe of positions and names (chr, pos, name)")
            return(0)
        }
    }
    clust = NULL
    if (overlayClusters) {
        if (is.null(sco@clust[1])) {
            print("ERROR: can't overlay clusters when clustering was not done on the input data")
            return(0)
        }
        else {
            clust = sco@clust
        }
    }
    num.rows <- length(cnToPlot) + 1
    if (add.legend) {
        num.rows <- num.rows + 1
    }
    textScale = 1
    axisPosScale = 1
    if (biggerText) {
        textScale = 1.4
        axisPosScale = 0.9
    }
    spacing = 1
    scale = 1
    if (num.rows == 2) {
        spacing = 1.5
        scale = 1.5
    }
    if (num.rows == 3) {
        spacing = 1
        scale = 1
    }
    height <- 8.5 * (num.rows/5)
    width <- 3.7
    if (!is.null(outputFile)) {
        pdf(file = outputFile, width = width, height = height, 
            bg = "white")
    }
    numClusters = 0
    if (!is.null(clust)) {
        numClusters = max(clust$cluster.assignments)
    }
    for (d in 1:dimensions) {
        name = sampleNames[d]
        vafs = getOneSampleVafs(vafs.merged, name, numClusters)
        par(mfcol = c(num.rows, 1), mar = c(0.5, 3/spacing, 1, 
            1.5/spacing), oma = c(3/spacing, 0.5, 4/spacing, 
            0.5), mgp = c(3, 1, 0))
        densities = densityData[[d]]$densities
        factors = densityData[[d]]$factors
        peakPos = densityData[[d]]$peakPos
        peakHeights = densityData[[d]]$peakHeights
        maxDepth = densityData[[d]]$maxDepth
        maxDensity = densityData[[d]]$maxDensity
        scalingFactor = 1/maxDensity
        plot.default(x = c(1:10), y = c(1:10), ylim = c(0, 1.1), 
            xlim = c(0, 100), axes = FALSE, ann = FALSE, col = "#00000000", 
            xaxs = "i", yaxs = "i")
        rect(0, 0, 100, 1.1, col = "#00000011", border = NA)
        axis(side = 2, at = c(0, 2), labels = c("", ""), las = 1, 
            cex.axis = 0.6 * textScale, hadj = 0.6 * textScale, 
            lwd = 0.5, lwd.ticks = 0.5, tck = -0.01)
        colors = c("#1C366099", "#67B32E99", "#F4981999", "#E5242099")
        density.curve.width <- 4
        for (i in cnToPlot) {
            if (!(is.null(densities[[i]])) & (!showHistogram | 
                (i != 2))) {
                lines(densities[[i]]$x, scalingFactor * factors[[i]], 
                  col = colors[i], lwd = density.curve.width)
                if (length(peakHeights[[i]]) > 0) {
                  ppos = c()
                  if (onlyLabelHighestPeak) {
                    ppos = which(peakHeights[[i]] == max(peakHeights[[i]]))
                  }
                  else {
                    ppos = which((peakHeights[[i]] == max(peakHeights[[i]])) & 
                      (peakHeights[[i]] > minimumLabelledPeakHeight))
                  }
                  if (length(ppos) > 0) {
                    text(x = peakPos[[i]][ppos], y = (scalingFactor * 
                      peakHeights[[i]][ppos]) + 1.7, labels = signif(peakPos[[i]][ppos], 
                      3), cex = 0.7, srt = 0, col = colors[[i]])
                  }
                }
            }
            else if (showHistogram & (i == 2)) {
                v = vafs[which(vafs$cleancn == 2 & vafs$adequateDepth == 
                  1), ]
                frequencies <- data.frame(x = v$vaf, row.names = NULL, 
                  stringsAsFactors = NULL)
                bin.width <- 2.5
                num.breaks <- ceiling(100/bin.width) + 1
                breaks <- unlist(lapply(0:(num.breaks - 1), function(x) 100 * 
                  x/(num.breaks - 1)))
                h <- hist(v$vaf, breaks = breaks, plot = FALSE)
                h$density <- h$density/max(h$density)
                plot(h, add = TRUE, freq = FALSE, col = "white", 
                  border = "black")
            }
        }
        model.style <- 4
        model.style <- 1
        model.width <- density.curve.width/2
        individual.model.style <- 3
        individual.model.width <- density.curve.width/2
        if (!(is.null(clust))) {
            maxFitDensity <- max(clust$fit.y[d, ])
            lines(clust$fit.x, clust$fit.y[d, ]/maxFitDensity, 
                type = "l", col = "grey50", lty = model.style, 
                lwd = model.width)
            if (overlayIndividualModels == TRUE) {
                for (i in 1:numClusters) {
                  lines(clust$fit.x, clust$individual.fits.y[[i]][d, 
                    ]/maxFitDensity, type = "l", col = "grey50", 
                    lty = individual.model.style, lwd = individual.model.width)
                }
            }
            if (highlightsOnHistogram) {
                if (!is.null(positionsToHighlight)) {
                  addpts = merge(vafs, positionsToHighlight, 
                    by.x = c("chr", "st"), by.y = c("chr", "st"))
                  for (i in 1:length(addpts$vaf)) {
                    if (addpts$name[i] != "") {
                      vaf <- addpts$vaf[i]
                      nearest.indx <- which(unlist(lapply(clust$fit.x, 
                        function(x) abs(x - vaf))) == min(abs(clust$fit.x - 
                        vaf)))[1]
                      vaf.y <- clust$fit.y[d, nearest.indx]/maxFitDensity
                      label <- as.character(addpts$name[i])
                      cex <- 1
                      text(x = vaf, y = vaf.y, label = "*", cex = cex)
                      text(x = vaf, y = vaf.y + 0.1, label = label, 
                        cex = cex)
                    }
                  }
                }
            }
        }
        lcol = colors[cnToPlot]
        lty = c(1, 1, 1, 1)
        lwd = c(2, 2, 2, 2)
        pchs = c(NA, NA, NA, NA)
        pt.bgs = lcol
        leg = c("1 Copy", "2 Copies", "3 Copies", "4 Copies")
        leg = leg[cnToPlot]
        if ((length(cnToPlot) == 1) & (cnToPlot[1] == 2)) {
            if (showHistogram == FALSE) {
                lty = c(1)
                lwd = c(2)
                pchs = c(NA)
                pt.bgs = lcol
            }
            else {
                lcol = "black"
                lty = c(0)
                lwd = c(0)
                pchs = c(22)
                pt.bgs = "white"
            }
        }
        if (!(is.null(clust))) {
            leg = c(leg, "Model Fit")
            lcol = c(lcol, "grey50")
            lty = c(lty, model.style)
            lwd = c(lwd, model.width)
            pt.bgs = c(pt.bgs, "grey50")
            pchs = c(pchs, NA)
            if (overlayIndividualModels == TRUE) {
                leg = c(leg, "Component Fits")
                lcol = c(lcol, "grey50")
                lty = c(lty, individual.model.style)
                lwd = c(lwd, 2)
                pt.bgs = c(pt.bgs, "grey50")
                pchs = c(pchs, NA)
            }
        }
        legend(x = "topright", lwd = lwd, lty = lty, legend = leg, 
            col = lcol, bty = "n", cex = ((0.6/scale) * textScale), 
            y.intersp = 1.25, pch = pchs, pt.bg = pt.bgs)
        axis(side = 3, at = c(0, 20, 40, 60, 80, 100), labels = c(0, 
            20, 40, 60, 80, 100), cex.axis = ((0.6/scale) * textScale), 
            lwd = 0.5, lwd.ticks = 0.5, padj = (((scale * 3.5) - 
                3.5) + 1.4) * (1/textScale), tck = -0.05)
        mtext("Variant Allele Frequency", adj = 0.5, padj = -3.1 * 
            (1/textScale), cex = 0.6 * (textScale * 0.8), side = 3)
        mtext("Density (a.u.)", side = 2, cex = 0.6 * (textScale * 
            0.8), padj = -4.2 * (axisPosScale))
        if (showTitle) {
            title = ""
            if (is.null(sampleNames[d])) {
                title = "Clonality Plot"
            }
            else {
                title = paste(sampleNames[d], "Clonality Plot", 
                  sep = " ")
            }
            mtext(title, adj = 0.5, padj = -5 * (1/textScale), 
                cex = 0.65 * textScale, side = 3)
        }
        if (showCopyNumberScatterPlots) {
            for (i in cnToPlot) {
                v = vafs[which(vafs$cleancn == i & vafs$adequateDepth == 
                  1), ]
                drawScatterPlot(v, highlightSexChrs, positionsToHighlight, 
                  colors, i, maxDepth, highlightsHaveNames, overlayClusters, 
                  scale, textScale, axisPosScale, labelOnPlot = FALSE, 
                  highlightCnPoints = highlightCnPoints)
                axis(side = 1, at = c(0, 20, 40, 60, 80, 100), 
                  labels = c(0, 20, 40, 60, 80, 100), cex.axis = (0.6/scale) * 
                    textScale, lwd = 0.5, lwd.ticks = 0.5, padj = ((-scale * 
                    5) + 5 - 1.4) * (1/textScale), tck = -0.05)
                if (length(cnToPlot) < 2 & highlightsHaveNames) {
                  addHighlightLegend(v, positionsToHighlight, 
                    scale)
                }
                else {
                  if (highlightsHaveNames) {
                    print("WARNING: highlighted point naming is only supported when plotting only CN2 regions (cnToPlot=c(2))")
                    print("Instead labeling directly on plot")
                  }
                }
            }
        }
    }
    if (!is.null(outputFile)) {
        devoff <- dev.off()
    }
}

getOneSampleVafs <- function (vafs.merged, name, numClusters) {
    common = c("chr", "st", "adequateDepth")
    a = vafs.merged[, common]
    prefix = paste("^", name, ".", sep = "")
    cols = grepl(prefix, names(vafs.merged))
    header = sub(prefix, "", names(vafs.merged)[cols])
    vafs = cbind(a, vafs.merged[, cols])
    names(vafs) = c(common, header)
    if (numClusters > 0) {
        vafs = cbind(vafs, vafs.merged$cluster)
        names(vafs)[length(vafs)] = "cluster"
    }
    return(vafs)
}

drawScatterPlot <- function (data, highlightSexChrs, positionsToHighlight, colors, 
    cn, maxDepth, highlightsHaveNames, overlayClusters, scale = 1, 
    textScale = 1, axisPosScale = 1, labelOnPlot = FALSE, highlightCnPoints = FALSE) {
    cex.points = 1/scale
    ptcolor = colors[cn]
    circlecolor = substr(colors[cn], 1, 7)
    plot.default(x = -10000, y = 1, log = "y", type = "p", pch = 19, 
        cex = 0.4, col = "#00000000", xlim = c(0, 100), ylim = c(5, 
            maxDepth * 3), axes = FALSE, ann = FALSE, xaxs = "i", 
        yaxs = "i")
    addPoints <- function(data, color, highlightSexChrs, pch = 16, 
        cex = 0.75, highlightCnPoints) {
        outlineCol = rgb(0, 0, 0, 0.1)
        data.sex = NULL
        data.cn = NULL
        if (highlightSexChrs) {
            data.sex = data[(data$chr == "X" | data$chr == "Y"), 
                ]
            data = data[!(data$chr == "X" | data$chr == "Y"), 
                ]
        }
        if (highlightCnPoints) {
            data.cn = data[data$chr == "CN", ]
            data = data[!(data$chr == "CN"), ]
        }
        points(data$vaf, data$depth, type = "p", pch = pch, cex = cex, 
            col = color)
        if (!(is.null(data.sex))) {
            points(data.sex$vaf, data.sex$depth, type = "p", 
                pch = pch + 1, cex = cex, col = color)
        }
        if (!(is.null(data.cn))) {
            points(data.cn$vaf, data.cn$depth, type = "p", pch = 15, 
                cex = cex + 0.5, col = "black")
            points(data.cn$vaf, data.cn$depth, type = "p", pch = 15, 
                cex = cex, col = "yellow")
        }
    }
    if (length(data[, 1]) > 0) {
        if (any(grepl(pattern = "^cluster$", names(data))) & 
            overlayClusters & cn == 2) {
            numClusters = max(data[, c("cluster")], na.rm = T)
            cols = getClusterColors(numClusters)
            p = data[data$cluster == 0, ]
            if (length(p[, 1]) > 0) {
                addPoints(p, rgb(0, 0, 0, 0.2), highlightSexChrs, 
                  cex = (cex.points * 0.6), highlightCnPoints = highlightCnPoints)
            }
            for (i in 1:numClusters) {
                p = data[data$cluster == i, ]
                addPoints(p, cols[i], highlightSexChrs, cex = cex.points, 
                  highlightCnPoints = highlightCnPoints)
            }
        }
        else {
            addPoints(data[data$chr != "CN", ], ptcolor, highlightSexChrs, 
                cex = cex.points, highlightCnPoints = highlightCnPoints)
            addPoints(data[data$chr == "CN", ], "black", highlightSexChrs, 
                pch = 15, cex = cex.points + 0.5, highlightCnPoints = highlightCnPoints)
            addPoints(data[data$chr == "CN", ], "yellow", highlightSexChrs, 
                pch = 15, cex = cex.points, highlightCnPoints = highlightCnPoints)
        }
        if (!(is.null(positionsToHighlight))) {
            addpts = merge(data, positionsToHighlight, by.x = c("chr", 
                "st"), by.y = c("chr", "st"))
            xs <- c()
            ys <- c()
            labels <- c()
            if (length(addpts[, 1]) > 0) {
                if (highlightsHaveNames) {
                  if (labelOnPlot) {
                    xs <- addpts$vaf[addpts$name != ""]
                    ys <- addpts$depth[addpts$name != ""]
                    labels <- addpts$name[addpts$name != ""]
                    num.labels <- length(labels)
                    suppressPackageStartupMessages(library(TeachingDemos))
                    text(x = addpts$vaf, y = addpts$depth, label = rep("*", 
                      length(addpts$vaf)))
                    if (num.labels > 0) {
                      xs <- spread.labs(xs, mindiff = 4, min = xs)
                      ys <- spread.labs(ys + 50, mindiff = 4, 
                        min = ys)
                      text(x = xs, y = ys, label = labels)
                    }
                  }
                  else {
                    for (i in 1:length(addpts$vaf)) {
                      if (addpts$name[i] != "") {
                        text(addpts$vaf[i], addpts$depth[i], 
                          labels = i, cex = 0.5)
                      }
                      else {
                        text(addpts$vaf[i], addpts$depth[i], 
                          labels = "*")
                      }
                    }
                  }
                }
                else {
                  addPoints(addpts, col = "#555555FF", highlightSexChrs, 
                    cex = cex.points, highlightCnPoints = highlightCnPoints)
                }
            }
        }
    }
    axis(side = 2, las = 1, tck = 0, lwd = 0, cex.axis = 1.2/(scale * 
        2) * textScale, hadj = 0.5/scale)
    for (i in 2:length(axTicks(2) - 1)) {
        lines(c(-1, 101), c(axTicks(2)[i], axTicks(2)[i]), col = "#00000022")
    }
    rect(-1, 5, 101, axTicks(2)[length(axTicks(2))] * 1.05, col = "#00000011", 
        border = NA)
    cnpos = axTicks(2)[length(axTicks(2)) - 1]
    points(x = c(95), y = cnpos, type = "p", pch = 19, cex = 3/scale, 
        col = circlecolor)
    text(c(95), y = cnpos, labels = c(cn), cex = 1/scale, col = "#FFFFFFFF")
    mtext("Tumor Coverage", side = 2, cex = 0.6 * (textScale * 
        0.8), padj = -4.2 * (axisPosScale))
}

getClusterColors <- function (numClusters) {
    suppressPackageStartupMessages(library(RColorBrewer))
    cols = brewer.pal(8, "Dark2")
    if (numClusters > 8) {
        cols = c(cols, brewer.pal(9, "Set1"))
    }
    if (numClusters > 20) {
        cols = rep(cols, ceiling(length(cols)/20))
    }
    cols = cols[1:numClusters]
    for (i in 1:length(cols)) {
        z = col2rgb(cols[i])
        cols[i] = rgb(z[1], z[2], z[3], 150, maxColorValue = 255)
    }
    return(cols)
}
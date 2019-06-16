##################################
##
## helper functions to render data for the Riboswitches pipeline report
##
##################################
suppressPackageStartupMessages({
  library(MASS)
  library(htmltools)
  library(knitr)
  library(DT)
  library(rjson)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(formattable)
  library(viridis)
  library(RColorBrewer)
  library(motifStack)
  library(MotIV)
  library(reshape2)
  library(igraph)
  library(WriteXLS)
})

#' loadGlobalVars: read configuration from bpipe vars, and assign them to the
#'  global environment.
#'
#' @param f - a file defining multiple variables for reporting to run.
#'
#' @return A set of variables that are mentioned in input file, e.g.
#'         PROJECT="/data/projects/projectID"
#'         QC="/data/projects/projectID/qc"
#'         FASTQ_PROCESSOR_RIBOSW_OUT="/data/projects/projectID/results/ribosw"
#'
loadGlobalVars <- function(f="report.conf") {
  sapply(readLines(f), function(x) {
    x <- unlist(strsplit(x, "=", fixed=T))
    assign(x[1], x[2], envir=.GlobalEnv)
  })

  invisible(NULL)
}

#' Read the config file and extract the patterns
#' Read and parse the config file containing the patterns to search for, and return
#' them in a named vector that can be used from R.
#'
#' @param f input config file where to read the json data from
#' @return named vector of PCRE
getPatternsFromConfig <- function(f="ribosw.json") {
  x <- fromJSON(paste(readLines(f), collapse=""))
  close(file(f))
  sapply(x$Patterns, function(x) {
    y <- x$RE
    names(y) <- x$Name
    y
  })
}

#' Read the config file and extract the targets
#' Read and parse the config file containing the sample description, and return
#' them in a data.frame that can be used from R.
#'
#' @param f input config file where to read the json data from
#' @return named vector of PCRE
getTargetsFromConfig <- function(f="ribosw.json") {
  x <- fromJSON(paste(readLines(f), collapse=""))
  close(file(f))
  targets <- as.data.frame(t(
    sapply(x$Targets, function(x) {
      c(projectID  = x$projectID,
        sampleID   = x$sampleID,
        sampleName = x$sampleName,
        group      = x$group,
        treatment  = x$treatment,
        subject    = x$subject,
        filePrefix = x$filePrefix)
    })
  ))

  # do some verifications there're no duplicated entries
  if(any(duplicated(targets$sampleID))   ||
     any(duplicated(targets$sampleName)) ||
     any(duplicated(targets$filePrefix)))
    stop("Duplicated sample ids/names/files in the Targets section of ribosw.json\n",
         "Will stop rendering the report. Check config file try rendering again.")

  # return targets
  targets
}

#' Read the config file and extract the contrasts
#' Read and parse the config file containing the contrasts to do the DE, and return
#' them in a data.frame that can be used from R.
#'
#' @param f input config file where to read the json data from
#' @return named vector of PCRE
getContrastsFromConfig <- function(f="ribosw.json") {
  x <- fromJSON(paste(readLines(f), collapse=""))
  close(file(f))
  as.data.frame(t(
    sapply(x$Contrasts, function(x) {
      c(Name   = x$Name,
        groupA = x$groupA,
        groupB = x$groupB,
        x      = paste0(x$Name, "=", x$groupA, "-", x$groupB))
    })
  ))
}

#' Read the config file and extract the model matrix
#' Read and parse the config file and return the model matrix to use with edgeR.
#'
#' @param f input config file where to read the json data from
#' @return named vector of PCRE
getMMatrixFromConfig <- function(f="ribosw.json") {
  x <- fromJSON(paste(readLines(f), collapse=""))
  close(file(f))
  x$ModelMatrix
}

#' Read the config file and extract the spike in to normalize against
#' Read and parse the config file and return the Spike in to use to normalize
#' counts.
#'
#' @param f input config file where to read the json data from
#' @return Spikein name
getSpikeinFromConfig <- function(f="ribosw.json") {
  x <- fromJSON(paste(readLines(f), collapse=""))
  close(file(f))
  x$UseSpikein
}

#' Get density of points in 2 dimensions.
#' From Kamil Slowikowski (http://slowkow.com/notes/ggplot2-color-by-density/)
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @return The density within each square.
#' @examples
#'   dat$density <- get.density(dat$x, dat$y)
#'   ggplot(dat) + geom_point(aes(x, y, color=density)) + scale_color_viridis()
get.density <- function(x, y, n=100) {
  ok <- is.finite(x) & is.finite(y)
  if(any(ok)) {
    dens <- try(MASS::kde2d(x=x[ok], y=y[ok], n=n), silent=TRUE)
    if(class(dens) == "try-error")
      dens <- try(MASS::kde2d(x=x[ok], y=y[ok], n=n, h=10), silent=TRUE)  # try to provide a bandwidth
    if(class(dens) == "try-error")
      return(0)
    ix <- findInterval(x, dens$x)  # may fail if x, y are not finite
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  } else {
    return(0) # no density
  }
}

#' panel.smooth for a pairs() plot
#' Modified version of the base function panel.smooth, in order to accomodate a
#' density color palette.
#' @param x vector of numbers corresponding to the x coordinates
#' @param y vector of numbers corresponding to the y coordinates
#' @param ... rest of the parameters that go to points() and lines()
panel.smooth.dens <- function (x, y, bg=NA, pch=par("pch"), cex=1,
                               col.smooth="red", span=2/3, iter=3, ...)
{
  ok <- is.finite(x) & is.finite(y)
  if(any(ok) && sum(ok) > 1) {
    dcols <- densCols(x=x[ok], y=y[ok], colramp=viridis, nbin=100)
    points(x[ok], y[ok], pch=pch, col=dcols, bg=bg, cex=cex)
    lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
  }
}

#' panel.hist for a pairs() plot
#' Taken literally from examples(pairs)
#' @param x vector of numbers
panel.hist <- function(x, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
}

#' helper to calculate the nucleotide (ACTG) frequencies per position
#' @param x character vector with ACTG nucleotides
#' @return matrix with ACTG counts per column
countsMat <- function(x) {

  # split
  x.split <- do.call(rbind, strsplit(sub("_", "", x), ""))
  rownames(x.split) <- x

  # and count
  apply(x.split, 2, function(x) {
    c(A=sum(x == "A"),
      C=sum(x == "C"),
      G=sum(x == "G"),
      T=sum(x == "T"))
  })
}

#' plot network from an edger's results table
# @param x a character vector of sequences
# @return a matrix with the hamming distances between sequences
hamming <- function(x) {
  sapply(x, function(a)
    sapply(x, function(b)
      sum(unlist(strsplit(a, "")) != unlist(strsplit(b, "")))
    )
  )
}

#' plot network from an edger's results table
# @param x a vector of patterns
# @param fc a vector of fc
# @param ... extra parms passed to plot()
plot.net.edger.results <- function(x, fc, legend.pos="topright", ...) {

  if(length(x) < 2) return()

  d <- hamming(x)
  am <- ifelse(d == 1, 1, 0) # remove distant sequences and autoreferences

  # make graph 
  net <- graph_from_adjacency_matrix(am, mode="undirected")

  # calculate the node color based on the FC
  pal <- rev(brewer.pal(11, "Spectral"))
  fc <- fc[match(names(V(net)), x)]
  cuts <- cut(fc, breaks=seq(from=-max(abs(fc)), to=max(abs(fc)), length.out=11), include.lowest=TRUE, right=FALSE)
  V(net)$color <- pal[cuts]
  V(net)$size  <- 6

  # detect communities and plot
  if(all(am == 0)) {    # there are no edges between vertex
    plot(net, layout=layout_with_fr, vertex.label.cex=.5,  # col=V(net)$color,
         ...)
  } else {
    ceb <- cluster_edge_betweenness(net) 
    ceb.communities <- communities(ceb)
    plot(ceb, net, layout=layout_with_fr, vertex.label.cex=.5, col=V(net)$color,
         mark.groups=ceb.communities[sapply(ceb.communities, length) > 1],
         mark.col="#00000000", mark.border="black",
         ...)
  }
  legend(legend.pos, fill=pal, legend=levels(cuts), bty="n", title="logFC", cex=2/3)
}

#' plot aligned motifs from an edger's results table
# @param pfms a list of pfms (see motifStack)
# @param fc a vector of fc
plot.piled.motifs <- function(pfms, fc, legend.pos="bottomleft") {

  # calculate the node color based on the FC
  pal <- rev(brewer.pal(11, "Spectral"))
  cuts <- cut(fc, breaks=seq(from=-max(abs(fc)), to=max(abs(fc)), length.out=11), include.lowest=TRUE, right=FALSE)
  col.fc <- pal[cuts]

  # convert the pfm matrix to a phylog object
  pfmList2matrixList <- function(pfms) {
      m <- lapply(pfms, pfm2pwm)
      names(m) <- unlist(lapply(pfms, function(.ele) .ele@name))
      m
  }

  # cluster and align motifs
  jaspar.scores <- MotIV::readDBScores(file.path(find.package("MotIV"), "extdata", "jaspar2010_PCC_SWU.scores"))
  d <- motifDistances(pfmList2matrixList(pfms), DBscores=jaspar.scores)
  hc <- motifHclust(d, method="average")
  fms <- pfms[hc$order]
  pfms <- DNAmotifAlignment(pfms)
  phylog <- hclust2phylog(hc)

  # plot
  motifPiles(phylog, pfms, r.anno=0.02, col.anno=list(fc=col.fc))
  legend(legend.pos, fill=pal, legend=levels(cuts), bty="n", title="logFC", cex=2/3)
}

#' Compute PI-value based on FC and p-value
#' PI-value is a score suggested in REF, which combines log2FC and FDR to assess significance.
#' It basically transforms the pvalue (-log10) by penalizing or enhancing it based on the log2FC.
#' Check the publication and supplemetary materials for a more detailed explanation.
#' REF: https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btr671
# @param pval a vector of p-values
# @param log2fc a vector of log2 tranformed fc
PIscore <- function(pval, log2fc) {
  list(score=log2fc * -log10(pval),   # the score
       pval =pval ^ log2fc)           # the transformed pval derived from the new score
}


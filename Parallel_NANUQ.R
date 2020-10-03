Parallel_quartetTable_byTree<-function (trees, taxonnames = NULL, 
                                          epsilon = 0, random = 0,
                                          parallel=FALSE,RAM_Gigs=1.0) {
  if (random < 0) 
    stop("Parameter 'random' must be non-negative.")
  random = round(random)
  if (is.null(taxonnames)) {
    taxonnames = trees[[1]]$tip.label
    message("Using taxa that appear on first gene tree.")
  }
  taxonnames = sort(taxonnames)
  nt = length(trees)
  N = length(taxonnames)
  if (random > 0) 
    M = random
  else M = choose(N, 4)
  Q = matrix(0, M, N + 4)
  qnames = c("12|34", "13|24", "14|23", "1234")
  colnames(Q) = c(taxonnames, qnames)
  warnMissing = 0
  if (random == 0) {
    m = 0
    for (i in 1:(N - 3)) {
      for (j in (i + 1):(N - 2)) {
        for (k in (j + 1):(N - 1)) {
          for (l in (k + 1):N) {
            m = m + 1
            Q[m, c(i, j, k, l)] = 1
          }
        }
      }
    }
  }
  else {
    i = 1
    while (i <= random) {
      q = sample(N, size = 4, replace = FALSE)
      row = integer(N)
      row[q] = 1
      j = 1
      while (j < i) {
        if (identical(Q[j, 1:N], row)) 
          j = i + 1
        else j = j + 1
      }
      if (j == i) {
        Q[i, 1:N] = row
        i = i + 1
      }
    }
  }
  
  
  
  message("Counting occurrences of displayed quartets for ", 
          M, " four-taxon subsets of ", N, " taxa across ", nt, 
          " gene trees.")
  start_time <- Sys.time()
  if(parallel==TRUE){
    options(future.globals.maxSize = RAM_Gigs * 1024^3) #1024^3 = ca. 1G RAM
    fopts <- options()[c("future.globals.maxSize")]
    array.out<-future_sapply(1:nt,FUN=function(numt,Q) {
      ## Set 'fopts' options (here fopts is a global variable)
      options(fopts)
      ## Show that those options are set
      cat(sprintf("test = %s\n", getOption("test")))
      t = trees[[numt]]
      if (is.null(t$edge.length)) {
        t = compute.brlen(t, 1)
      }
      else {
        if (sum(t$edge.length < 0) > 0) 
          stop("Error: Negative branch length in tree")
      }
      zeros = which(t$edge.length <= epsilon)
      t$edge.length[] = 1
      t$edge.length[zeros] = 0
      d = cophenetic.phylo(t)
      for (m in 1:M) {
        tax = as.character(taxonnames[which(Q[m, 1:N] == 
                                              1)])
        if (all(tax %in% colnames(d))) {
          a = d[tax[1], tax[2]] + d[tax[3], tax[4]]
          b = d[tax[1], tax[3]] + d[tax[2], tax[4]]
          c = d[tax[1], tax[4]] + d[tax[2], tax[3]]
          z = sort(c(a, b, c))
          if (z[1] == z[2]) {
            Q[m, "1234"] = Q[m, "1234"] + 1
          }
          else {
            if (z[1] == a) {
              Q[m, "12|34"] = Q[m, "12|34"] + 1
            }
            else {
              if (z[1] == b) {
                Q[m, "13|24"] = Q[m, "13|24"] + 1
              }
              else {
                Q[m, "14|23"] = Q[m, "14|23"] + 1
              }
            }
          }
        }
        else warnMissing = 1
      }
      return(Q[,c("12|34", "13|24", "14|23", "1234")])
      # return(Q)
    },Q=Q,simplify = "array")
  }else{
    array.out<-sapply(1:nt,FUN=function(numt,Q) {
      t = trees[[numt]]
      if (is.null(t$edge.length)) {
        t = compute.brlen(t, 1)
      }
      else {
        if (sum(t$edge.length < 0) > 0) 
          stop("Error: Negative branch length in tree")
      }
      zeros = which(t$edge.length <= epsilon)
      t$edge.length[] = 1
      t$edge.length[zeros] = 0
      d = cophenetic.phylo(t)
      for (m in 1:M) {
        tax = as.character(taxonnames[which(Q[m, 1:N] == 
                                              1)])
        if (all(tax %in% colnames(d))) {
          a = d[tax[1], tax[2]] + d[tax[3], tax[4]]
          b = d[tax[1], tax[3]] + d[tax[2], tax[4]]
          c = d[tax[1], tax[4]] + d[tax[2], tax[3]]
          z = sort(c(a, b, c))
          if (z[1] == z[2]) {
            Q[m, "1234"] = Q[m, "1234"] + 1
          }
          else {
            if (z[1] == a) {
              Q[m, "12|34"] = Q[m, "12|34"] + 1
            }
            else {
              if (z[1] == b) {
                Q[m, "13|24"] = Q[m, "13|24"] + 1
              }
              else {
                Q[m, "14|23"] = Q[m, "14|23"] + 1
              }
            }
          }
        }
        else warnMissing = 1
      }
      return(Q[,c("12|34", "13|24", "14|23", "1234")])
      # return(Q)
    },Q=Q,simplify = "array")
  }
  
  array.summed<-rowSums(array.out, dims = 2)
  Q<-cbind(Q[,-which(names(as.data.frame(Q))%in%qnames)],array.summed)
  # Q<-array.out
  
  if (warnMissing == 1) 
    warning("Some taxa missing from some trees.")
  if ((N > 4) && (sum(rowSums(Q[, qnames, drop = FALSE]) == 
                      0) > 0)) 
    warning("Some 4-taxon set not present on any tree.")
  if (sum(Q[, "1234"]) > 0) 
    warning("Some quartets unresolved.")
  current_time = Sys.time()
  elapsedTime = as.numeric(difftime(current_time, start_time, 
                                    units = "mins"))
  # if (elapsedTime > 15) {
  message("Time to process quartets on gene trees was ", 
          elapsedTime, " minutes.")
  # }
  return(Q)
}

### implement parallel quartetTable in NANUQ

Parallel_NANUQ<-function (genedata, outfile = "NANUQdist", alpha = 0.05, beta = 0.95, 
                            taxanames = NULL, plot = TRUE, RAM_Gigs = 1.0) 
{
  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("Critical values alpha and beta must be numeric.")
  }
  if ("matrix" %in% class(genedata)) {
    pTable = genedata
    if (!is.null(taxanames)) {
      message("Ignoring argument taxanames since genedata supplied as quartet table.")
    }
  }
  else {
    if ("multiPhylo" %in% class(genedata)) {
      genetrees = genedata
    }
    else {
      if ("character" %in% class(genedata)) {
        genetrees <- read.tree(genedata)
        message(paste("Read", length(genetrees), "gene trees from file."))
      }
      else {
        stop("Data must be supplied as an object of type multiPhylo, character, or matrix.")
      }
    }
    if (is.null(taxanames)) {
      taxanames = genetrees[[1]]$tip.label
    }
    taxanames = sort(taxanames)
    if (length(taxanames) <= 25) {
      namelist = paste0(taxanames, collapse = ", ")
    }
    else {
      namelist = paste0(paste0(taxanames[1:25], collapse = ", "), 
                        ",...(see output table for full list)")
    }
    message("Analyzing ", length(taxanames), " taxa: ", 
            namelist)
    pTable = Parallel_quartetTable_byTree(genetrees, taxonnames=taxanames,
                                            parallel=TRUE,RAM_Gigs = RAM_Gigs)
    pTable = quartetTableResolved(pTable, omit = FALSE)
  }
  if (!("p_T3" %in% colnames(pTable))) {
    pTable = quartetTreeTestInd(pTable, model = "T3", lambda = 0)
  }
  if (!("p_star" %in% colnames(pTable))) {
    pTable = quartetStarTestInd(pTable)
  }
  NANUQdist(pTable, outfile, alpha, beta, plot)
  invisible(pTable)
}

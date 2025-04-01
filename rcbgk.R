# representing causal background knowledge

library(pcalg)
library(RBGL)

# ======== transform cbgk =========
criticalSet <- function(g, xlb, zlb){
  
  # given a chordal graph g, g may be disconnected
  # find the critical set of x w.r.t z
  # g is a graphNEL object
  # xlb is a vertex (label name), zlb is an array (label names)
  
  x <- which(g@nodes == xlb)
  z <- which(g@nodes %in% zlb)
  n <- length(g@nodes)
  S <- list(c(x, 0, 0))
  C <- matrix(0, 1, n)
  flag <- 1
  dlp <- 0
  while (length(S) != 0){
    dlp <- dlp +1
    e <- S[[1]]
    S[[1]] <- NULL
    if (sum(z == e[1]) == 0 && e[2] != 0){
      for (alpha in setdiff(g@edgeL[[e[1]]][[1]],e[2])){
        if (sum(g@edgeL[[e[2]]][[1]]==alpha)==0){
          if (e[2] == x){
            S[[length(S)+1]] <- c(alpha, e[1], e[1])
          }else{
            S[[length(S)+1]] <- c(alpha, e[1], e[3])
          }
        }
      }
    }else if (sum(z == e[1])==0){
      for (alpha in g@edgeL[[e[1]]][[1]]){
        S[[length(S)+1]] <- c(alpha, e[1], alpha)
      }
    }else {
      C[flag] <- e[3]
      flag <- flag + 1
    }
    # if chordless cycle presents, then the while loop will never end
    # DEAD loop prevent
    if (dlp > (n^2+1)) break
  }
  C <- unique(C[C!=0])
  return(g@nodes[C])
}

find.ancestors.and.chaincomp <- function(amat, x){
  # amat is a adj matrix of a cpdag
  # x is a given node, label name
  # this function attempts to find the set of ancestors, chain components, 
  # and critical set of each node in the chain component
  
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  labelNames <- colnames(amat)
  
  w.an <- c(x)
  res.an <- c()
  while (length(w.an) > 0){
    node <- w.an[1]
    w.an <- w.an[-1]
    an.tmp <- labelNames[amatDir[, node] == 1]
    w.an <- append(w.an, setdiff(setdiff(an.tmp, res.an), w.an))
    res.an <- append(res.an, node)
  }
  
  w.cp <- c(x)
  res.cp <- c()
  while (length(w.cp) > 0){
    node <- w.cp[1]
    w.cp <- w.cp[-1]
    cp.tmp <- labelNames[amatUdir[, node] == 1]
    w.cp <- append(w.cp, setdiff(setdiff(cp.tmp, res.cp), w.cp))
    res.cp <- append(res.cp, node)
  }
  
  return(list(an = res.an, cp = res.cp))
  
}

transCbgk <- function(cpdag, bgk){
  # cpdag, graphNEL obj, a cpdag
  # bgk, a matrix with 3 columns, (tail, head, type), could be integers or characters
  # integers are positions, characters are label names
  # where type == 1 means causal, type == 0 means non-causal
  # output, a list vectors, direct causal clauses
  # vertices are indexed by label names
  
  amat <- as(cpdag, 'matrix')
  n <- nrow(amat)
  labelNames <- colnames(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  ug <- as(amatUdir, 'graphNEL')
  
  res <- lapply(1:n, function(.) c())
  names(res) <- labelNames
  
  ancestor <- lapply(1:n, function(.) c())
  names(ancestor) <- labelNames
  
  cp <- lapply(1:n, function(.) c())
  names(cp) <- labelNames
  
  if (nrow(bgk) != 0){
    for (i in 1:nrow(bgk)){
      if (is.numeric(bgk[i, 1])) {x <- labelNames[bgk[i, 1]]} else {x <- bgk[i, 1]}
      if (is.numeric(bgk[i, 2])) {y <- labelNames[bgk[i, 2]]} else {y <- bgk[i, 2]}     
      t <- as.integer(bgk[i, 3])
      if (amatSkel[x, y] == 1){
        # x is adjacent to y
        if (t == 0) {res[[y]] <- append(res[[y]], list(c(x)))} else {res[[x]] <- append(res[[x]], list(c(y)))}

      }else{
        # x is not adjacent to y
        # find the critical set Cxy of x w.r.t. y
        cset <- criticalSet(cpdag, x, y)
        
        if (length(cset) != 0){
          if (t == 0){
            for (element in cset) res[[element]] <- append(res[[element]], list(c(x)))
          }else{
            res[[x]] <- append(res[[x]], list(cset))
          }
        }else{
          if (t == 1) res[[x]] <- append(res[[x]], list(c())) # contradiction
        }
      }
    }
  }
  return(Filter(Negate(is.null), res))
}

# ======== random generation ========
genValidDcc <- function(cpdag, num){
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  
  tails <- sample(labelNames, num, replace = TRUE)
  uniqTails <- table(tails)
  dcc <- lapply(1:length(uniqTails), function(.) c())
  names(dcc) <- names(uniqTails)
  
  for (i in 1:length(uniqTails)){
    x <- names(dcc)[i]
    adjx <- labelNames[amatSkel[, x] == 1]
    pax <- labelNames[amatDir[, x] == 1]
    chx <- labelNames[amatDir[x, ] == 1]
    
    if (all(adjx %in% pax)){
      next
    }
    xnum <- uniqTails[x]
    for (j in 1:xnum){
      ind <- sample(c(0,1), length(adjx), replace = TRUE)
      select <- adjx[ind == 1]
      newSelect <- c()
      if (all(select %in% pax)){
        can <- setdiff(adjx, pax)
        while (length(newSelect) == 0){
          newSelect <- can[sample(c(0,1), length(can), replace = TRUE) == 1]
        }
      }
      if (!all(select %in% chx)){
        if (runif(1) < 0.8) select <- setdiff(select, chx)
      }
      D <- union(select, newSelect)
      dcc[[i]] <- append(dcc[[i]], list(c(D)))
    }
    dcc[[i]] <- unique(Filter(length, dcc[[i]]))
  }
  dcc <- Filter(length, dcc)
  return(dcc)
}

genInvalidDcc <- function(cpdag, num, prob = 0.9){
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  
  tails <- sample(colnames(amat), num, replace = TRUE)
  uniqTails <- table(tails)
  dcc <- lapply(1:length(uniqTails), function(.) c())
  names(dcc) <- names(uniqTails)
  
  for (i in 1:length(uniqTails)){
    x <- names(dcc)[i]
    adjx <- labelNames[amatSkel[, x] == 1]
    pax <- labelNames[amatDir[, x] == 1]
    chx <- labelNames[amatDir[x, ] == 1]
    
    xnum <- uniqTails[x]
    for (j in 1:xnum){
      ind <- sample(c(0,1), length(adjx), replace = TRUE)
      select <- adjx[ind == 1]
      newSelect <- c()
      if (runif(1) < prob) newSelect <- sample(colnames(amat), 1)
      D <- union(select, newSelect)
      dcc[[i]] <- append(dcc[[i]], list(c(D)))
    }
    dcc[[i]] <- unique(Filter(length, dcc[[i]]))
  }
  dcc <- Filter(length, dcc)
  return(dcc)
}

genConsistentCbgk <- function(cpdag, num, type = 1){
  # type: 1, direct causal relations 
  #       2, non-causal relations
  #       3, causal relations
  
  amat.cpdag <- as(cpdag, 'matrix')
  amat.cpdag[amat.cpdag != 0] <- 1 
  res <- pcalg::pdag2dag(as(amat.cpdag, 'graphNEL'))
  if (res$success){
    dag <- res$graph
    amatDagT <- t(as(dag, 'matrix'))
    nn = nrow(amatDagT)
    labelNames <- colnames(amatDagT)
  }else{
    stop('the input cpdag is not extenable.')
  }
  
  if (type == 1){
    edgeList <- which(as(dag, 'matrix') != 0, arr.ind = TRUE)
    positions <- sample(nrow(edgeList), num)
    xx <- edgeList[positions, 1]
    yy <- edgeList[positions, 2]
    cbgk <- cbind(xx, yy, 1)
  }else if (type == 2){
    xx <- c()
    yy <- c()
    for (i in 1:nn){
      de <- possDe(amatDagT, i, type = "dag")
      non.ancestor <- setdiff(1:nn, de)
      yy <- c(yy, non.ancestor)
      xx <- c(xx, rep(i, length(non.ancestor)))
    }
    positions <- sample(length(xx), num)
    cbgk <- cbind(xx[positions], yy[positions], 0)
  }else if (type == 3){
    xx <- c()
    yy <- c()
    tt <- c()
    for (i in 1:nn){
      de <- possDe(amatDagT, i, type = "dag")
      yy <- c(yy, setdiff(de, i))
      xx <- c(xx, rep(i, length(de)-1 ))
      tt <- c(tt, rep(1, length(de)-1 ))
      
    }
    positions <- sample(length(xx), num)
    cbgk <- cbind(xx[positions], yy[positions], tt[positions])
  }else{
    stop('the input type is not implemented')
  }
  cbgk[, 1:2] <- labelNames[cbgk[, 1:2]]
  return(cbgk)
}

# ======== validity and consistency ========
isValid <- function(cpdag, dcc){
  # dcc: label names
  
  amat <- as(cpdag, 'matrix')
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  labelNames <- colnames(amat)

  if (length(dcc) != 0){
    for (i in 1:length(dcc)){
      x <- names(dcc)[i]
      for (j in 1:length(dcc[[i]])){
        D <- dcc[[i]][[j]]
        if ((length(D) == 0) || any(amatSkel[x, D] == 0) || all(amatDir[D, x] == 1)){
          return(FALSE)
        }
      }
    }
    
  }
  return(TRUE)
}

isConsistent <- function(cpdag, dcc){
  # dcc: label names
  
  amat <- as(cpdag, 'matrix')
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  labelNames <- colnames(amat)
  
  if (!isValid(cpdag, dcc)) 
    return(list(indicator = FALSE, PEO = NULL))
  
  indicator <- TRUE
  vset <- labelNames
  peo <- c()
  while (length(vset) > 0){
    pln <- restriction(cpdag, vset, dcc)$potentialLeaf
    if (length(pln) == 0){
      indicator <- FALSE
      peo <- NULL
      break
    }else{
      vset <- setdiff(vset, pln[1])
      peo <- c(peo, pln[1])
    }
  }
  return(list(indicator = indicator, PEO = peo))
}

# ======== construct MPDAG ========
constructMOC <- function(cpdag, x, dcc, simplified=TRUE){
  # dcc: label names
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  
  if (!isValid(cpdag, dcc)) 
    stop("the direct causal clauses are valid with respect to the CPDAG.")
  
  vset <- labelNames
  while (length(vset) > 0){
    pln <- restriction(cpdag, vset, dcc)$potentialLeaf
    # print(pln)
    if (length(pln) == 0){
      stop("the direct causal clauses are inconsistent with respect to the CPDAG.")
    }else if (length(pln) == 1 && pln == x){
      mocSet <- vset
      break
    }else{
      vset <- setdiff(vset, setdiff(pln, x)[1])
    }
  }
  adj = intersect(labelNames[amatUdir[, x] == 1], vset)
  if (simplified) return(adj)
  return(list(MOC = vset,  adj = adj))
}

constructMPDAG <- function(cpdag, dcc){
  # dcc: label names
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  new.pa <- lapply(labelNames, constructMOC, cpdag=cpdag, dcc=dcc)
  names(new.pa) <- labelNames
  new.pa <- Filter(length, new.pa)
  if (length(new.pa) != 0){
    for (i in 1:length(new.pa)){
      x <- names(new.pa)[i]
      p <- new.pa[[i]]
      amat[x, p] <- 0
    }
  }
  return(as(amat, 'graphNEL'))
}

constructMPDAG.fast <- function(cpdag, dcc){
  # only consider non-singleton chain components
  # and seperately consider each connected chain components
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  ug <- as(amatUdir, 'graphNEL')
  cc <- connectedComp(ug)
  for (flag in 1:length(cc)){
    vset <- cc[[flag]]
    rst <- restriction(cpdag, vset, dcc)
    new.pa <- lapply(vset, constructMOC, cpdag=rst$UdirSubGraph, dcc=rst$restricted)
    names(new.pa) <- vset
    new.pa <- Filter(length, new.pa)
    if (length(new.pa) != 0){
      for (i in 1:length(new.pa)){
        x <- names(new.pa)[i]
        p <- new.pa[[i]]
        amat[x, p] <- 0
      }
    }
  }
  return(as(amat, 'graphNEL'))
}

# ======== restriction ========
restriction <- function(cpdag, vset, dcc){
  # dcc: label names
  # return the undirected induced subgraph of cpdag over vset
  # also return the restriction of dcc on that subgraph
  
  if (length(vset) == 0)
    stop("vset should be non-empty.")
  
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  amatUdirSub <- as.matrix(amatUdir[vset, vset])
  colnames(amatUdirSub) <- rownames(amatUdirSub) <- vset
  usubg <- as(amatUdirSub, 'graphNEL')
  
  if (!isValid(cpdag, dcc)) 
    stop("the direct causal clauses are invalid with respect to the CPDAG.")

  restrictDcc <- dcc[names(dcc) %in% vset]
  if (length(restrictDcc) > 0){
    for (i in 1:length(restrictDcc)){
      x <- names(restrictDcc)[i]
      for (j in 1:length(restrictDcc[[i]])){
        D <- restrictDcc[[i]][[j]]
        
        # restrict x \tor D to G^*_u
        if (any(amatDir[x, D] == 1)){
          restrictDcc[[i]][j] <- list(c())
        }else{
          restrictDcc[[i]][[j]] <- intersect(restrictDcc[[i]][[j]], labelNames[amatUdir[x, ] == 1])
        }
        
        if (!all(restrictDcc[[i]][[j]] %in% vset))
          restrictDcc[[i]][j] <- list(c())
      }
      restrictDcc[[i]] <- unique(Filter(Negate(is.null), restrictDcc[[i]]))
    }
  }
  restrictDcc <- Filter(length, restrictDcc)
  
  # find potential leaf nodes
  candidates <- setdiff(colnames(amatUdirSub), names(restrictDcc))
  pln <- c()
  if (length(candidates) != 0){
    pln <- candidates[sapply(candidates, FUN = isComplete, amat = amatUdirSub)]
  }
  
  
  return(list(UdirSubGraph = usubg, restricted = restrictDcc, potentialLeaf = pln))
}

isComplete <- function(amat, node){
  # node: label name
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  adjNodesAmat <- as.matrix(amatSkel[which(amatSkel[node,] == 1), which(amatSkel[node,] == 1)])
  if (sum(adjNodesAmat == 0) == ncol(adjNodesAmat)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# ======== IDA ========
bgkIDA <- function(x, y, mcov, cpdag, mpdag = NULL, bgk = NULL, verbose = FALSE){
  
  # x, y can be integers (positions) or characters (label names)
  
  labelNames <- cpdag@nodes

  if (is.numeric(x)){
    xlb <- labelNames[x]
  }else{
    xlb <- x
    if (!(xlb %in% labelNames)) stop('x is not in the cpdag')
    x <- which(labelNames == xlb)
  } 
  
  if (is.numeric(y)){
    ylb <- labelNames[y]
  }else{
    ylb <- y
    if (!(ylb %in% labelNames)) stop('y is not in the cpdag')
    y <- which(labelNames == ylb)
  } 
  
  if (length(bgk) == 0){
    beta.hat <- ida(x, y, mcov, cpdag, method = "local", type = 'cpdag')
  }else{
    # determine whether the input cpdag is valid
    amat.cpdag <- ad.g <- wgtMatrix(cpdag)
    amat.cpdag[which(amat.cpdag != 0)] <- 1
    if (!isValidGraph(amat = amat.cpdag, type = 'cpdag')) {
      message("The input graph is not a valid cpdag. See function isValidGraph() for details.\n")
    }
    
    # trans cbgk to dcc
    dcc <- transCbgk(cpdag, bgk)
    if (length(mpdag) == 0){
      mpdag <- constructMPDAG(cpdag, dcc)
    }
    
    # determine whether the constructed mpdag is valid
    amat <- ad.g <- wgtMatrix(mpdag)
    amat[which(amat != 0)] <- 1
    if (!isValidGraph(amat = amat, type = 'pdag')) {
      message("The input graph is not a valid pdag. See function isValidGraph() for details.\n")
    }
    
    nl <- colnames(amat)
    stopifnot(!is.null(nl))
    amatSkel <- amat + t(amat)
    amatSkel[amatSkel != 0] <- 1
    
    # local method
    wgt.est <- (ad.g != 0)       # the adjacency matrix of the mpdag, transpose
    tmp <- wgt.est - t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp            # directed subgraph, transpose
    pa1 <- which(wgt.unique[x, ] != 0)       # definite parents
    if (y %in% pa1) {
      beta.hat <- 0
    }else{
      wgt.ambig <- wgt.est - wgt.unique      # undirected subgraph
      pa2 <- which(wgt.ambig[x, ] != 0)      # siblings
      if (verbose) 
        cat("\n\nx=", x, "y=", y, "\npa1=", pa1, "\npa2=", pa2, "\n")
      if (length(pa2) == 0){
        # calculate causal effect
        beta.hat <- lm.cov(mcov, y, c(x, pa1))
        if (verbose) 
          cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", beta.hat, "\n")
      }else{
        beta.hat <- NA
        ii <- 0
        pa2.f <- pa2           # out from x  
        pa2.t <- NULL          # into x
        if (isLocallyConsistent(cpdag, dcc, x, pa2.t, pa2.f)) {
          ii <- ii + 1
          beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1))
          if (verbose) 
            cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", beta.hat[ii], "\n")
        }
        for (i2 in seq_along(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (isLocallyConsistent(cpdag, dcc, x, pa2.t, pa2.f)){
            ii <- ii + 1
            if (y %in% pa2.t) {
              beta.hat[ii] <- 0
            }else {
              beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1, pa2[i2]))
              if (verbose) 
                cat("Fit - y:", y, "x:", c(x, pa1, pa2[i2]), "|b.hat=", beta.hat[ii], "\n")
            }
          }
        }
        if (length(pa2) > 1){
          for (i in 2:length(pa2)){
            pa.tmp <- combn(pa2, i, simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))){
              pa2.t <- pa.tmp[, j]
              pa2.f <- setdiff(pa2, pa2.t)
              if (isLocallyConsistent(cpdag, dcc, x, pa2.t, pa2.f)) {
                ii <- ii + 1
                if (y %in% pa2.t) {
                  beta.hat[ii] <- 0
                }else{
                  beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1, pa2.t))
                  if (verbose) 
                    cat("Fit - y:", y, "x:", c(x, pa1, pa2.t), "|b.hat=", beta.hat[ii], "\n")
                }
              }
            }
          }
        } 
      }
    }
  }
  unname(beta.hat)
}

isLocallyConsistent <- function(cpdag, dcc, x, pa2.t, pa2.f){
  # dcc: label names
  # x, pa2.t, pa2.f: integers
  
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  vset <- union(labelNames[amatUdir[, x] == 1], labelNames[x])
  pa2.t <- labelNames[pa2.t]
  pa2.f <- labelNames[pa2.f]
  x <- labelNames[x]
  
  dcc.new <- dcc
  if (length(pa2.t) != 0){
    for (p in pa2.t){
      dcc.new[[p]] <- append(dcc.new[[p]], list(c(x)))
    }
  }
  if (length(pa2.f) != 0){
    for (c in pa2.f){
      dcc.new[[x]] <- append(dcc.new[[x]], list(c(c)))
    }
  }
  
  rdcc <- restriction(cpdag, vset, dcc.new)
  indicator <- isConsistent(rdcc$UdirSubGraph, rdcc$restricted)
  return(indicator$indicator)
}

lm.cov <- function (C, y, x) {
  # borrowed from R packages pcalg
  # see, https://github.com/cran/pcalg
  
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}

# ======== from bgIDA ========

criticalSetOld <- function(g, x, z){
  # given a chordal graph g, g may be disconnected
  # find the critical set of x w.r.t z
  # g is a graphNEL object
  # x is a number, z is an array
  
  n <- length(g@nodes)
  S <- list(c(x, 0, 0))
  C <- matrix(0, 1, n)
  flag <- 1
  dlp <- 0
  while (length(S) != 0){
    dlp <- dlp +1
    e <- S[[1]]
    S[[1]] <- NULL
    if (sum(z == e[1])==0 && e[2] != 0){
      for (alpha in setdiff(g@edgeL[[e[1]]][[1]],e[2])){
        if (sum(g@edgeL[[e[2]]][[1]]==alpha)==0){
          if (e[2] == x){
            S[[length(S)+1]] <- c(alpha, e[1], e[1])
          }else{
            S[[length(S)+1]] <- c(alpha, e[1], e[3])
          }
        }
      }
    }else if (sum(z == e[1])==0){
      for (alpha in g@edgeL[[e[1]]][[1]]){
        S[[length(S)+1]] <- c(alpha, e[1], alpha)
      }
    }else {
      C[flag] <- e[3]
      flag <- flag + 1
    }
    # if chordless cycle presents, then the while loop wil never end
    # DEAD loop prevent
    if (dlp > (n^2+1)) break
  }
  C <- unique(C[C!=0])
  return(C)
}

find.ancestors.and.chaincomp.old <- function(amat, x){
  # amat is a adj matrix of a cpdag
  # x is a given node
  # this function attempts to find the set of ancestors 
  # and chain components in the chain component
  
  n <- nrow(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  
  w.an <- c(x)
  res.an <- c()
  while (length(w.an) > 0){
    node <- w.an[1]
    w.an <- w.an[-1]
    an.tmp <- which(amatDir[, node] == 1)
    w.an <- append(w.an, setdiff(setdiff(an.tmp, res.an), w.an))
    res.an <- append(res.an, node)
  }
  
  w.cp <- c(x)
  res.cp <- c()
  while (length(w.cp) > 0){
    node <- w.cp[1]
    w.cp <- w.cp[-1]
    cp.tmp <- which(amatUdir[, node] == 1)
    w.cp <- append(w.cp, setdiff(setdiff(cp.tmp, res.cp), w.cp))
    res.cp <- append(res.cp, node)
  }
  
  return(list(an = res.an, cp = res.cp))
  
}

add.bg.withlabel <- function(cpdag, xx = c(), yy = c()){
  # another implementation of add.bg
  # if the nodes in cpdag have labels
  # and xx, yy are the labels of nodes instead of the indices of nodes
  # please use this function instead of add.bg
  
  amat <- as(cpdag, 'matrix')
  n <- nrow(amat)
  labels <- colnames(amat)
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  ug <- as(amatUdir, 'graphNEL')
  # res.x and res.y store characters, i.e. nodes labels
  res.x <- c()
  res.y <- c()
  ancestor <- lapply(1:n, function(.) c())
  cp <- lapply(1:n, function(.) c())
  
  if (length(xx) != 0){
    for (i in 1:length(xx)){
      # x and y are labels
      x <- as.character(xx[i])
      y <- as.character(yy[i])
      pos.x <- which(labels == x)
      pos.y <- which(labels == y)
      if (amatSkel[x, y] == 1){
        # x is adjacent to y
        res.x <- c(res.x, x)
        res.y <- c(res.y, y)
      }else{
        # x is not adjacent to y
        # y is not an ancestor of x
        # for each z \in an(x, cpdag) \cap cp(y), y is not an ancestor of z
        # find the critical set of y w.r.t z
        
        if (is.null(ancestor[[pos.x]])){
          tmp <- find.ancestors.and.chaincomp.old(amat, pos.x)
          ancestor[[pos.x]] <- tmp$an
          cp[[pos.x]] <- tmp$cp
        }
        if (is.null(cp[[pos.y]])){
          tmp <- find.ancestors.and.chaincomp.old(amat, pos.y)
          ancestor[[pos.y]] <- tmp$an
          cp[[pos.y]] <- tmp$cp
        }
        zset <- intersect(ancestor[[pos.x]], cp[[pos.y]])
        if (length(zset) != 0){
          c <- criticalSetOld(ug, pos.y, zset)
          res.x <- c(res.x, labels[c])
          res.y <- c(res.y, rep(y, length(c)))
        }
      }
    }
  }
  pdag <- addBgKnowledge(cpdag, x = res.x, y = res.y)
}

dida <- function(x, y, mcov, graphEst, verbose = FALSE) {
  
  type = 'pdag'
  
  # graphEst is a maximal PDAG
  
  method <- "local"
  y.notparent = FALSE
  
  # graphEst <- addBgKnowledge(cpdag, x = bn[, 1], y = bn[, 2])
  amat <- ad.g <- wgtMatrix(graphEst)
  amat[which(amat != 0)] <- 1
  if (!isValidGraph(amat = amat, type = type)) {
    message("The input graph is not a valid ", type, ". See function isValidGraph() for details.\n")
  }
  nl <- colnames(amat)
  stopifnot(!is.null(nl))
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1
  if (method == "local") {
    wgt.est <- (ad.g != 0) # cpdag, t
    tmp <- wgt.est - t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp # directed subgraph, t
    pa1 <- which(wgt.unique[x, ] != 0) # definite pa
    if (y %in% pa1) {
      beta.hat <- 0
    }
    else {
      wgt.ambig <- wgt.est - wgt.unique  # undirected subgraph
      pa2 <- which(wgt.ambig[x, ] != 0)  # sib
      if (verbose) 
        cat("\n\nx=", x, "y=", y, "\npa1=", pa1, "\npa2=", 
            pa2, "\n")
      if (length(pa2) == 0) {
        # calculate causal effect
        beta.hat <- lm.cov(mcov, y, c(x, pa1))
        if (verbose) 
          cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", 
              beta.hat, "\n")
      }
      else {
        beta.hat <- NA
        ii <- 0
        pa2.f <- pa2
        pa2.t <- NULL
        if (hasExtensionNew(amat, amatSkel, x, pa1, pa2.t, 
                            pa2.f)) {
          ii <- ii + 1
          beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1))
          if (verbose) 
            cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", 
                beta.hat[ii], "\n")
        }
        for (i2 in seq_along(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (hasExtensionNew(amat, amatSkel, x, pa1, pa2.t, 
                              pa2.f)) {
            ii <- ii + 1
            if (y %in% pa2.t) {
              beta.hat[ii] <- 0
            }
            else {
              beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1, 
                                                pa2[i2]))
              if (verbose) 
                cat("Fit - y:", y, "x:", c(x, pa1, pa2[i2]), 
                    "|b.hat=", beta.hat[ii], "\n")
            }
          }
        }
        if (length(pa2) > 1) 
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2, i, simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))) {
              pa2.t <- pa.tmp[, j]
              pa2.f <- setdiff(pa2, pa2.t)
              if (hasExtensionNew(amat, amatSkel, x, pa1, 
                                  pa2.t, pa2.f)) {
                ii <- ii + 1
                if (y %in% pa2.t) {
                  beta.hat[ii] <- 0
                }
                else {
                  beta.hat[ii] <- lm.cov(mcov, y, c(x, 
                                                    pa1, pa2.t))
                  if (verbose) 
                    cat("Fit - y:", y, "x:", c(x, pa1, 
                                               pa2.t), "|b.hat=", beta.hat[ii], 
                        "\n")
                }
              }
            }
          }
      }
    }
  }
  unname(beta.hat)
}

hasExtensionNew <- function(amat, amatSkel, x, pa1, pa2.t, pa2.f){
  # borrowed from R packages pcalg
  # see, https://github.com/cran/pcalg
  
  res <- !has.new.coll.or.cycle(amat, amatSkel, x, pa1, pa2.t, pa2.f)
  res
}

has.new.coll.or.cycle <- function(amat,amatSkel, x, pa1, pa2.t, pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## or directed triangle containing x
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if a new collider or a directed triangle is introduced 
  
  res <- FALSE
  
  # check v-structure
  if (length(pa2.t) > 0 && !all(is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if not, there is a new collider
    if (length(pa1) > 0 && !all(is.na(pa1))) {
      res <- min(amatSkel[pa1, pa2.t]) == 0 ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if (!res && length(pa2.t) > 1) {
      A2 <- amatSkel[pa2.t,pa2.t]
      diag(A2) <- 1
      res <- min(A2) == 0 ## TRUE if new collider
    }
  }
  if (!res && length(pa2.f) > 0 && !all(is.na(pa2.f))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    A <- amat-t(amat)
    A[A < 0] <- 0
    ## find parents of pa2.f
    cA <- colSums(A[pa2.f,,drop = FALSE])
    papa <- setdiff(which(cA != 0), x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa) > 0)
      res <- min(amatSkel[x,papa]) == 0 ## TRUE if new collider
  }
  
  ## checking direcrted triangle containing X
  if (!res ) {
    ## consider here pa = pa1 U pa2.t, cd = pa2.f U cd(amat)
    ## cd should not point towards pa
    
    ## adding new orientations to directed subgraph
    A <- amat-t(amat)
    A[A < 0] <- 0
    nA <- A
    nA[x, pa2.t] <- 1
    nA[pa2.t, x] <- 0
    nA[x, pa2.f] <- 0
    nA[pa2.f, x] <- 1
    
    # check whether cd point towards pa in nA
    cd <- which(nA[, x] != 0)
    pa <- which(nA[x, ] != 0) 
    if (length(cd) > 0 && length(pa) > 0){
      subA <- nA[pa, cd]
      res <- max(subA) == 1 ## TRUE if exists cd->pa
    }
  }
  res
}

# ======== utility functions ========

bgkIDA.semilocal <- function(x, y, mcov, cpdag, mpdag = NULL, bgk = NULL, verbose = FALSE){
  # x, y can be integers (positions) or characters (label names)
  
  labelNames <- cpdag@nodes
  
  if (is.numeric(x)){
    xlb <- labelNames[x]
  }else{
    xlb <- x
    if (!(xlb %in% labelNames)) stop('x is not in the cpdag')
    x <- which(labelNames == xlb)
  } 
  
  if (is.numeric(y)){
    ylb <- labelNames[y]
  }else{
    ylb <- y
    if (!(ylb %in% labelNames)) stop('y is not in the cpdag')
    y <- which(labelNames == ylb)
  } 
  
  if (length(bgk) == 0){
    beta.hat <- ida(x, y, mcov, cpdag, method = "local", type = 'cpdag')
  }else{
    # determine whether the input cpdag is valid
    amat.cpdag <- ad.g <- wgtMatrix(cpdag)
    amat.cpdag[which(amat.cpdag != 0)] <- 1
    if (!isValidGraph(amat = amat.cpdag, type = 'cpdag')) {
      message("The input graph is not a valid cpdag. See function isValidGraph() for details.\n")
    }
    
    # trans cbgk to dcc
    dcc <- transCbgk(cpdag, bgk)
    if (length(mpdag) == 0){
      mpdag <- constructMPDAG(cpdag, dcc)
    }
    
    # determine whether the constructed mpdag is valid
    amat <- ad.g <- wgtMatrix(mpdag)
    amat[which(amat != 0)] <- 1
    if (!isValidGraph(amat = amat, type = 'pdag')) {
      message("The input graph is not a valid pdag. See function isValidGraph() for details.\n")
    }
    
    nl <- colnames(amat)
    stopifnot(!is.null(nl))
    amatSkel <- amat + t(amat)
    amatSkel[amatSkel != 0] <- 1
    
    # local method
    wgt.est <- (ad.g != 0)       # the adjacency matrix of the mpdag, transpose
    tmp <- wgt.est - t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp            # directed subgraph, transpose
    pa1 <- which(wgt.unique[x, ] != 0)       # definite parents
    if (y %in% pa1) {
      beta.hat <- 0
    }else{
      wgt.ambig <- wgt.est - wgt.unique      # undirected subgraph
      pa2 <- which(wgt.ambig[x, ] != 0)      # siblings
      if (verbose) 
        cat("\n\nx=", x, "y=", y, "\npa1=", pa1, "\npa2=", pa2, "\n")
      if (length(pa2) == 0){
        # calculate causal effect
        beta.hat <- lm.cov(mcov, y, c(x, pa1))
        if (verbose) 
          cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", beta.hat, "\n")
      }else{
        beta.hat <- NA
        ii <- 0
        pa2.f <- pa2           # out from x  
        pa2.t <- NULL          # into x
        if (isLocallyConsistent(cpdag, dcc, x, pa2.t, pa2.f)) {
          ii <- ii + 1
          beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1))
          if (verbose) 
            cat("Fit - y:", y, "x:", c(x, pa1), "|b.hat=", beta.hat[ii], "\n")
        }
        for (i2 in seq_along(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (isLocallyConsistent(cpdag, dcc, x, pa2.t, pa2.f)){
            ii <- ii + 1
            if (y %in% pa2.t) {
              beta.hat[ii] <- 0
            }else {
              beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1, pa2[i2]))
              if (verbose) 
                cat("Fit - y:", y, "x:", c(x, pa1, pa2[i2]), "|b.hat=", beta.hat[ii], "\n")
            }
          }
        }
        if (length(pa2) > 1){
          for (i in 2:length(pa2)){
            pa.tmp <- combn(pa2, i, simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))){
              pa2.t <- pa.tmp[, j]
              pa2.f <- setdiff(pa2, pa2.t)
              if (isLocallyConsistent(cpdag, dcc, x, pa2.t, pa2.f)) {
                ii <- ii + 1
                if (y %in% pa2.t) {
                  beta.hat[ii] <- 0
                }else{
                  beta.hat[ii] <- lm.cov(mcov, y, c(x, pa1, pa2.t))
                  if (verbose) 
                    cat("Fit - y:", y, "x:", c(x, pa1, pa2.t), "|b.hat=", beta.hat[ii], "\n")
                }
              }
            }
          }
        } 
      }
    }
  }
  unname(beta.hat)
}

isGloballyConsistent <- function(cpdag, dcc, x, pa2.t, pa2.f){
  # dcc: label names
  # x, pa2.t, pa2.f: integers
  
  amat <- as(cpdag, 'matrix')
  labelNames <- colnames(amat)
  amatDir <- amat - t(amat)
  amatDir[which(amatDir < 0)] <- 0
  amatUdir <- amat - amatDir
  ug <- as(amatUdir, 'graphNEL')
  
  pa2.t <- labelNames[pa2.t]
  pa2.f <- labelNames[pa2.f]
  x <- labelNames[x]
  
  dcc.new <- dcc
  if (length(pa2.t) != 0){
    for (p in pa2.t){
      dcc.new[[p]] <- append(dcc.new[[p]], list(c(x)))
    }
  }
  if (length(pa2.f) != 0){
    for (c in pa2.f){
      dcc.new[[x]] <- append(dcc.new[[x]], list(c(c)))
    }
  }

  indicator <- isConsistent(ug, dcc)
  return(indicator$indicator)
}

## ======== evaluation ========
effect_distance <- function(x, y, tdag, effect_pair){
  # compute true causal effects pair
  true.mat <- as(tdag, 'matrix')
  mcov <- trueCov(tdag)
  amatDag <- true.mat
  amatDag[amatDag != 0] <- 1
  
  pa <- which(true.mat[, x] != 0)
  de <- possDe(t(amatDag), x, type = "dag")
  if (y %in% de){
    true.eff <- lm.cov(mcov, y, c(x, pa))
  }else{
    true.eff <- 0
  }
  
  
  # evaluate
  distmat <- as.matrix(dist(c(true.eff, effect_pair)))
  cmse <- sum(distmat[2:nrow(distmat), 1]^2)/(nrow(distmat) - 1)
  
  return(c(min(distmat[2:nrow(distmat), 1]), max(distmat[2:nrow(distmat), 1]), cmse) )
}

int_mse <- function(x, y, tdag, effect_pair){
  # compute true causal effects pair
  true.mat <- as(tdag, 'matrix')
  mcov <- trueCov(tdag)
  amatDag <- true.mat
  amatDag[amatDag != 0] <- 1
  
  pa <- which(true.mat[, x] != 0)
  de <- possDe(t(amatDag), x, type = "dag")
  if (y %in% de){
    true.eff <- lm.cov(mcov, y, c(x, pa))
  }else{
    true.eff <- 0
  }
  
  
  # evaluate
  ub <- max(effect_pair)
  lb <- min(effect_pair)
  if (true.eff > ub){
    return(true.eff - ub)
  }else if (true.eff < lb){
    return(lb - true.eff)
  }else{
    return(0)
  }
}


















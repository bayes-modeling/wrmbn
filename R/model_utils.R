
structure_learning <- function(data, number_layers, bl, wl, cluster, slearning_algo = "tabu", number_bootstrap = 100, debug = FALSE) {

  if(debug) {
    cat(file = stderr(), "Learning structure by algorithm", slearning_algo, ", number bootstrap", number_bootstrap, "\n")
  }

  data_kpi <- data
  #Start learning
  dyn.dag = NULL
  if(slearning_algo == "tabu") {
    dyn.dag = bnlearn::tabu(data_kpi, blacklist = bl, whitelist = wl, max.iter = 10000, optimized = FALSE)
  } else if(slearning_algo == "fast.iamb") {
    dyn.dag = bnlearn::fast.iamb(data_kpi, blacklist = bl, whitelist = wl, cluster = cluster)
  } else if(slearning_algo == "inter.iamb") {
    dyn.dag = bnlearn::fast.iamb(data_kpi, blacklist = bl, whitelist = wl, cluster = cluster)
  } else if(slearning_algo == "iamb") {
    dyn.dag = bnlearn::iamb(data_kpi, blacklist = bl, whitelist = wl, cluster = cluster)
  } else if(slearning_algo == "hc") {
    dyn.dag = bnlearn::hc(data_kpi, blacklist = bl, whitelist = wl, max.iter = 10000, optimized = FALSE)
  } else if(slearning_algo == "gs") {
    dyn.dag = bnlearn::gs(data_kpi, blacklist = bl, whitelist = wl, cluster = cluster)
  } else {
    return(NULL)
  }

  #Bootstrap
  if(debug) {
    cat(file = stderr(), "Start bootstrap structure by algorithm", slearning_algo, ", number bootstrap", number_bootstrap, "\n")
  }

  dyn.bootstrap = custom_bootstrap(data_kpi, number_bootstrap, slearning_algo, bl, wl, cluster, debug)
  threshold = attr(dyn.bootstrap, "threshold")

  #Pruning arcs
  dyn.avg = bnlearn::averaged.network(dyn.bootstrap, threshold = threshold)

  undirected_arcs <- bnlearn::undirected.arcs(dyn.avg)
  if(dim(undirected_arcs)[1] != 0) {
    for(i in 1:length(undirected_arcs[, 1])) {
      dyn.avg <- bnlearn::drop.arc(dyn.avg, from = undirected_arcs[i, 1], to = undirected_arcs[i, 2])
    }
  }

  #Learning parameters
  fitted <- bnlearn::bn.fit(dyn.avg, data, cluster = cluster)

  result <- list()
  result[["dyn"]] <- dyn.dag
  result[["dyn.avg"]] <- dyn.avg
  result[["fitted"]] <- fitted
  result[["bootstrap"]] <- dyn.bootstrap
  result[["threshold"]] <- threshold
  if(debug) {
    cat(file = stderr(), "Learning structure by algorithm", slearning_algo, ", number bootstrap", number_bootstrap, "completed", "\n")
  }
  return(result)
}

remove_unknown_arc <- function(arcs, variables) {
  remove_index <- c()
  for(i in 1:length(arcs[, 1])) {
    fr <- arcs[i, 1]
    to <- arcs[i, 2]

    if(length(which(fr %in% variables)) == 0 || length(which(to %in% variables)) == 0) {
      remove_index <- c(remove_index, i)
    }
  }
  if(length(remove_index) > 0) {
    arcs <- arcs[-remove_index, ]
  }

  if(length(arcs[, 1]) > 0) {
    rownames(arcs) <- 1:length(arcs[, 1])
    return(arcs)
  } else {
    return(NULL)
  }
}

custom_bootstrap <- function(data, R, algorithm, blacklist, whitelist, cluster, debug = FALSE) {
  if(debug) {
    cat(file = stderr(), paste("Boostrap with R = ", R, " and algorithm = ", algorithm), "\n")
  }

  dyn.str = bnlearn::boot.strength(data, cluster = cluster, R = R, algorithm = algorithm, algorithm.args = list(blacklist = blacklist, whitelist = whitelist))

  if(debug) {
    cat(file = stderr(), "Boostrap with R = ", R, " and algorithm = ", algorithm, " completed", "\n")
  }

  return(dyn.str)
}


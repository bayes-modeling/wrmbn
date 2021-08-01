structure_filter <- function(known_structure, desire_layers, dynamic_variables,
                             is_blacklist_internal, is_blacklist_other, is_variable_only,
                             rho, corr_threshold,
                             custom_blacklist, custom_whitelist) {


  whitelist <- filter_whitelist(known_structure, desire_layers, dynamic_variables,
                                is_blacklist_internal, is_variable_only, rho, corr_threshold,
                                custom_blacklist, custom_whitelist)

  blacklist <- filter_blacklist(desire_layers, dynamic_variables, is_blacklist_internal, is_blacklist_other,
                                custom_blacklist, whitelist)

  return(list(whitelist = whitelist, blacklist = blacklist))

}

filter_whitelist <- function(known_structure, desire_layers, dynamic_variables,
                             is_blacklist_internal, is_variable_only, rho, corr_threshold,
                             custom_blacklist, custom_whitelist) {
  graph_full <- NULL

  if(!is.null(known_structure)) {
    graph_full <- get_known_structure(desire_layers, dynamic_variables, known_structure, is_variable_only)
    graph_full <- filter_unknow_arc(graph_full, dynamic_variables)
  } else {
    graph_full <- generate_fully_connected_graph(dynamic_variables, desire_layers)
  }

  graph_full <- filter_by_correlation(graph_full, rho, corr_threshold)

  if(is_blacklist_internal) {
    graph_full <- filter_by_blacklist_internal(graph_full, dynamic_variables, desire_layers)
  }

  if(!is.null(custom_blacklist)) {
    graph_full <- filter_by_custom_blacklist(graph_full, custom_blacklist$datapath)
  }

  if(!is.null(custom_whitelist)) {
    graph_full <- add_custom_whitelist(graph_full, custom_whitelist$datapath)
  }

  return(graph_full)
}

filter_blacklist <- function(desire_layers, dynamic_variables, is_blacklist_internal, is_blacklist_other,
                             custom_blacklist, whitelist = NULL) {


  if(is_blacklist_other && !is.null(whitelist)) {
    return(blacklist_all(dynamic_variables, whitelist, desire_layers))
  } else {
    blacklist <- get_reverse_time_arcs(desire_layers, dynamic_variables)

    if(is_blacklist_internal) {
      blacklist_internal <- get_all_internal_arcs(desire_layers, dynamic_variables)
      blacklist <- rbind(blacklist, blacklist_internal)
    }

    if(!is.null(custom_blacklist)) {
      blacklist_custom <- get_custom_blacklist(custom_blacklist$datapath)
      blacklist <- rbind(blacklist, blacklist_custom)
    }

    blacklist <- filter_unknow_arc(blacklist, dynamic_variables)

    return(blacklist)
  }


}

get_all_discrete_values <- function(data, discrete_variables) {
  return(unique(data[, discrete_variables]))
}

normalize_data <- function(data, type = "mean_normalization", custom_normalizers = NULL) {
  cols <- colnames(data)

  normalizers <- list()

  for(col in cols) {
    values <- data[, col]

    if(is.numeric(values)) {
      values <- values[which(!is.na(values))]
      normalizer <- NULL

      if(is.null(custom_normalizers)
         | is.null(custom_normalizers[[col]])) {
        normalizer <- get_normalizer(values)
      } else {
        normalizer <- custom_normalizers[[col]]
      }

      normalizers[[col]] <- normalizer
      if(type == "mean_normalization") {
        data[, col] <- (data[, col] - normalizer[["mean"]]) / (normalizer[["max"]] - normalizer[["min"]])
      } else if(type == "min_max") {
        data[, col] <- (data[, col] - normalizer[["min"]]) / (normalizer[["max"]] - normalizer[["min"]])
      } else if(type == "standardisation") {
        data[, col] <- (data[, col] - normalizer[["mean"]]) / sd[["min"]]
      }
    }
  }

  result <- list()
  result[["data"]] <- data
  result[["normalizers"]] <- normalizers

  return(result)
}

get_normalizer <- function(values) {
  normalizer <- list()

  normalizer[["mean"]] <- mean(values)
  normalizer[["min"]] <- min(values)
  normalizer[["sd"]] <- sd(values)
  normalizer[["max"]] <- max(values)

  return(normalizer)
}

extract_cross_validation_loss <- function(algorithm.cv) {
  loss <- c()
  for(i in 1:length(algorithm.cv)) {
    run <- algorithm.cv[[i]]

    for(j in 1:length(run)) {
      fold <- run[[j]]
      loss <- c(loss, fold$loss)
    }
  }

  return(loss)
}

group_levels <- function(data, data_levels) {
  cols <- colnames(data)
  group_level <- matrix(0, nrow = length(data_levels), ncol = length(cols))
  colnames(group_level) <- cols
  rownames(group_level) <- data_levels

  for(i in 1:length(cols)) {
    value_levels <- unique(data[, i])
    for(level in value_levels){
      if(level %in% data_levels) {
        count <- sum(data[, i] == level)
        proportion <- round(count / nrow(data), 3) * 100
        group_level[level, i] = proportion
      }
    }
  }

  return(group_level)
}

group_continuous <- function(data) {
  group_result <- matrix(0, nrow = ncol(data), ncol = 5)
  colnames(group_result) <- c("Lower 95", "Lower 85", "Mean", "Upper 85", "Upper 95")
  rownames(group_result) <- colnames(data)
  for(i in 1:ncol(data)) {
    values <- data[, i]
    error_95 <- round(qt(0.975,df=length(values)-1)*sd(values)/sqrt(length(values)), 5)
    error_85 <- round(qt(0.925,df=length(values)-1)*sd(values)/sqrt(length(values)), 5)

    value_mean <- round(mean(values), 5)
    lower_85 <- value_mean - error_85
    upper_85 <- value_mean + error_85

    lower_95 <- value_mean - error_95
    upper_95 <- value_mean + error_95

    group_result[i, 1] <- lower_95
    group_result[i, 2] <- lower_85
    group_result[i, 3] <- value_mean
    group_result[i, 4] <- upper_85
    group_result[i, 5] <- upper_95
  }

  return(group_result)
}

remove_bias_variables <- function(data, pattern, threshold) {
  data_length = dim(data)[1]

  remove_col_index <- c()
  for(col in colnames(data)) {
    number_zeros = sum(as.character(data[, col]) == pattern)
    if((number_zeros / data_length) >= threshold) {
      remove_col_index <- c(remove_col_index, which(col == colnames(data)))
    }
  }
  if(length(remove_col_index) != 0) {
    return(data[, -remove_col_index])
  }
  return(data)
}

convert_data_type <- function(data, continuous_columns, discrete_columns) {
  data_cols <- colnames(data)
  for(continuous_col in continuous_columns) {
    time_cols <- grep(continuous_col, data_cols, value = TRUE)
    for(time_col in time_cols) {
      data[, time_col] <- as.numeric(data[, time_col])
    }
  }

  return(data)
}

convert_variables_to_factor <- function(data) {
  data_cols <- colnames(data)
  data_factors <- list()
  for(data_col in data_cols) {
    if(is.character(data[, data_col])) {
      data_factor <- as.factor(data[, data_col])
      data_factors[[data_col]] <- levels(data_factor)
      data[, data_col] <- as.numeric(data_factor)
    } else if(is.integer(data[, data_col])) {
      data_factor <- as.factor(data[, data_col])
      data_factors[[data_col]] <- levels(data_factor)
      data[, data_col] <- as.numeric(data_factor)
    } else if(is.factor(data[, data_col])){
      data_factors[[data_col]] <- levels(data[, data_col])
      data[, data_col] <- as.numeric(data[, data_col])
    }
  }

  return(list(data = data, data_factors = data_factors))
}

convert_to_numeric <- function(discrete_data) {
  col_names <- colnames(discrete_data)

  for(col_name in col_names) {
    col_data <- as.double(as.factor(discrete_data[, col_name]))
    discrete_data[, col_name] <- col_data
  }

  return(discrete_data)
}

get_dynamic_variables <- function(variables, layers) {
  dynamic_variables <- c()

  for(variabel in variables) {
    dynamic_variables <- c(dynamic_variables, paste(variabel, 1:layers, sep = "_"))
  }

  return(dynamic_variables)
}

init_blacklist <- function(total_variables, number_of_layer) {
  bl <- data.frame(from = character(), to = character())
  for(i in number_of_layer:1) {
    froms <- grep(paste(i, "$", sep = ""), total_variables, value = TRUE)
    if(i > 1) {
      for(j in 1:i-1) {
        tos <- grep(paste(j, "$", sep = ""), total_variables, value = TRUE)
        empty <- expand.grid(from = froms, to = tos, stringsAsFactors = FALSE)
        bl <- rbind(bl, empty)
      }
    } else {
      empty.t1 <- expand.grid(from = froms, to = froms, stringsAsFactors = FALSE)
      bl <- rbind(bl, empty.t1)
    }
  }
  return(bl)
}

deduplicate_discrete <- function(discrete_data, threshold) {
  discrete_data_numeric <- convert_to_numeric(discrete_data) + 0.01
  discrete_data_numeric <- dedup(discrete_data_numeric, threshold)

  return(discrete_data[, colnames(discrete_data_numeric)])
}

deduplicate_continuous <- function(discrete_continuous, threshold) {
  continuous_data_numeric <- discrete_continuous + 0.01
  continuous_data_numeric <- dedup(continuous_data_numeric, threshold)

  return(discrete_continuous[, colnames(continuous_data_numeric)])
}

generate_fully_connected_graph <- function(dynamic_variables, number_layers) {
  graph <- data.frame(from = character(), to = character())
  for(i in 1:(number_layers - 1)) {
    from_i = grep(paste0(i, "$"), dynamic_variables, value = TRUE)
    for(j in (i+1): number_layers) {
      to_j = grep(paste0(j, "$"), dynamic_variables, value = TRUE)
      graph_from_i <- expand.grid(from = from_i, to = to_j, stringsAsFactors = FALSE)
      graph <- rbind(graph, graph_from_i)
    }
  }

  colnames(graph) <- c("from", "to")
  return(graph)
}

filter_by_correlation <- function(graph, cor_matrix, threshold) {
  remove_arcs <- c()

  for(i in 1:nrow(graph)) {
    from <- graph[i, 1]
    to <- graph[i, 2]

    if(abs(cor_matrix[from, to]) <= threshold) {

      remove_arcs <- c(remove_arcs, i)
    }
  }
  if(length(remove_arcs) > 0) {
    graph <- graph[-remove_arcs, ]
  }

  return(graph)
}

filter_by_blacklist_internal <- function(graph, variables, number_layers) {
  for(i in 1:(number_layers - 1)) {
    layer_variables <- grep(paste0(i, "$"), variables, value = TRUE)
    internal_arcs <- generate_internal_arcs(layer_variables)
    graph <- get_graph_diff(graph, internal_arcs)

  }

  return(graph)
}

filter_by_custom_blacklist <- function(graph, blacklist_file_path) {
  custom_blacklist <- read.csv(blacklist_file_path)

  return(get_graph_diff(graph, custom_blacklist))
}

get_custom_blacklist <- function(blacklist_file_path) {
  return(read.csv(blacklist_file_path))
}

get_all_internal_arcs <- function(number_layers, variables) {
  internal_arcs <- data.frame(from = character(), to = character())
  for(i in 1:(number_layers - 1)) {
    internal_layer_variables <- grep(paste0(i, "$"), variables, value = TRUE)

    internal_layer_arcs <- generate_internal_arcs(internal_layer_variables)
    internal_arcs <- rbind(internal_arcs, internal_layer_arcs)
  }

  return(internal_arcs)
}

blacklist_all <- function(variables, whitelist, number_layers) {
  all_arcs <- expand.grid(from = variables, to = variables, stringsAsFactors = FALSE)
  blacklist <- get_graph_diff(all_arcs, whitelist)

  return(blacklist)
}

get_reverse_time_arcs <- function(number_layers, variables) {
  bl <- data.frame(from = character(), to = character())
  for(i in number_layers:1) {
    froms <- grep(paste(i, "$", sep = ""), variables, value = TRUE)
    if(i > 1) {
      for(j in 1:i-1) {
        tos <- grep(paste(j, "$", sep = ""), variables, value = TRUE)
        empty <- expand.grid(from = froms, to = tos, stringsAsFactors = FALSE)
        bl <- rbind(bl, empty)
      }
    }
  }
  return(bl)
}

generate_internal_arcs <- function(layer_variables) {
  layer_graph <- expand.grid(from = layer_variables, to = layer_variables, stringsAsFactors = FALSE)
  return(layer_graph)
}

get_known_structure <- function(number_layers, all_variables, structure, is_variable_only) {
  known_structure <- data.frame(from = character(), to = character())

  if(!is.null(structure) && nrow(structure) > 0) {
    if(number_layers > 1 && is_variable_only) {
      froms <- structure$from
      tos <- structure$to

      for(i in 1:number_layers) {
        from_i <- paste(froms, i, sep = "_")
        if(i < number_layers) {
          #Generate structure between layers
          for(j in (i + 1):number_layers) {

            to_j <- paste(tos, j, sep = "_")
            known_i_to_j <- data.frame(from = from_i, to = to_j)
            known_structure <- rbind(known_structure, known_i_to_j)
          }
        } else {
          #Generate structure within last layer
          to_i <- paste(tos, i, sep = "_")
          known_i_to_i <- data.frame(from = from_i, to = to_i)
          known_structure <- rbind(known_structure, known_i_to_i)
        }

      }
      self_structure <- init_self_structure(all_variables, number_layers)

      known_structure <- rbind(known_structure, self_structure)
    } else {
      known_structure <- structure
    }

  }

  return(known_structure)

}

get_variable_name <- function(variable) {
  variables <- strsplit(variable, "_")[[1]]
  return(variables[length(variables) - 1])
}

get_variable_layer <- function(variable) {
  variables <- strsplit(variable, "_")[[1]]
  return(as.numeric(variables[length(variables)]))
}


init_self_structure <- function(all_variables, number_layers) {
  wl <- c("from", "to")
  for(variable_1 in all_variables) {
    for(variable_2 in all_variables) {
      if(variable_1 != variable_2) {
        variable_1_name <- get_variable_name(variable_1)
        variable_2_name <- get_variable_name(variable_2)

        variable_1_layer <- get_variable_layer(variable_1)
        variable_2_layer <- get_variable_layer(variable_2)

        if(variable_1_name == variable_2_name
           && variable_1_layer < variable_2_layer) {
          wl <- rbind(wl, c(variable_1, variable_2))
        }
      }
    }
  }

  wl <- wl[-1, ]
  colnames(wl) <- c("from", "to")
  return(wl)
}


get_group <- function(variables) {
  groups <- c()
  for(variable in variables) {
    variable_elems <- strsplit(variable, "_")[[1]]
    group = variable_elems[length(variable_elems)]
    if(length(variable_elems) > 1) {
      groups <- c(groups, group)
    } else {
      groups <- c(groups, 0)
    }

  }

  return(groups)
}

filter_unknow_arc <- function(graph, all_variables) {
  remove_index <- c()
  for(i in 1:nrow(graph)) {
    from <- graph$from[i]
    to <- graph$to[i]

    if((length(which(from %in% all_variables)) == 0
        || length(which(to %in% all_variables)) == 0)
       || from == to) {
      remove_index <- c(remove_index, i)
    }
  }

  if(length(remove_index) > 0) {
    graph <- graph[-remove_index, ]
  }
  return(graph)
}

add_custom_whitelist <- function(graph, whitelist_file_path) {
  custom_whitelist <- read.csv(whitelist_file_path)
  graph <- rbind(graph, custom_whitelist)
  graph <- graph[-duplicated(graph), ]
  return(graph)
}

get_graph_diff <- function(graph1, graph2) {
  remove_index <- c()
  for(i in 1:nrow(graph1)) {
    from_1 <- graph1[i, "from"]
    to_1 <- graph1[i, "to"]
    for(j in 1:nrow(graph2)) {
      from_2 <- graph2[j, "from"]
      to_2 <- graph2[j, "to"]

      if(from_1 == from_2 && to_1 == to_2) {
        remove_index <- c(remove_index, i)
        break
      }
    }
  }

  if(length(remove_index) > 0) {
    return(graph1[-remove_index, ])
  }

  return(graph1)
}

generate_sector_equation <- function(sector_value, sector_variable, range) {

  equation <- paste("(", sector_variable, " >= ", sector_value - range, " & ", sector_variable, " <= ", sector_value + range, ")", sep = "")

  return(equation)
}

generate_dist_node <- function(nodes) {
  node_strings <- c()
  for(node in nodes) {
    node_strings <- c(node_strings, paste("\"", node, "\"", sep = ""))
  }
  equation <- paste("nodes = c(", paste(node_strings, collapse = ","), ")", sep = "")
  return(equation)
}

get_kpi_score <- function(kpi_value, kpi_quantile, score_table) {
  if(is.na(kpi_value)) return(-1)
  for(i in 2:length(kpi_quantile)) {
    if(kpi_value <= kpi_quantile[i]) {
      return(score_table[i])
    }
  }
  return(0)
}

get_anomally_kpi <- function(kpi_score, threshold) {
  kpi_values <- as.numeric(kpi_score[1, ])
  return(colnames(kpi_score)[which(-1 < kpi_values & kpi_values <= threshold)])
}



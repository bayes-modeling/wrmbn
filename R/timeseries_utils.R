# Title     : TODO
# Objective : TODO
# Created by: ADMIN
# Created on: 6/19/2021
reverse_normalized_data <- function(data, type = NULL, normalizers = NULL) {
  if(is.null(type) | is.null(normalizers)) return(data)

  cols <- colnames(data)

  for(col in cols) {
    values <- data[, col]

    if(is.numeric(values)) {
      normalizer <- NULL

      if(!is.null(normalizers[[col]])) {
        normalizer <- normalizers[[col]]

        if(type == "mean_normalization") {
          data[, col] <- data[, col] * (normalizer[["max"]] - normalizer[["min"]]) + normalizer[["mean"]]
        } else if(type == "min_max") {
          data[, col] <- data[, col] * (normalizer[["max"]] - normalizer[["min"]]) + normalizer[["min"]]
        } else if(type == "standardisation") {
          data[, col] <- data[, col] * sd[["min"]] + normalizer[["mean"]]
        }
      }
    }
  }

  return(data)
}


reverse_constructed_data <- function(data, desire_layers) {

  all_variables <- colnames(data)
  variable_without_time <- unique(gsub("*_[1-9]$", "", all_variables))
  data_matrix <- matrix(0, nrow = nrow(data) + desire_layers - 1, ncol = length(variable_without_time))
  data_matrix <- data.frame(data_matrix)
  colnames(data_matrix) <- variable_without_time

  if(nrow(data) == 1) {
    for(variable in variable_without_time) {
      variable_times <- grep(variable, all_variables, value = TRUE)
      data_matrix[, variable] <- as.numeric(data[1, variable_times])
    }
  } else {
    for(variable in variable_without_time) {
      variable_times <- grep(variable, all_variables, value = TRUE)
      variable_1 <- data[, variable_times[1]]
      variable_last <- data[, variable_times[desire_layers]]
      variable_values <- c(variable_1,  variable_last[(length(variable_last) - desire_layers + 2):length(variable_last)])
      data_matrix[, variable] <- variable_values
    }
  }

  return(data_matrix)
}
reconstruct_timeseries_data <- function(data, number_layers, na_omit = TRUE) {
  cols <- colnames(data)

  new_cols <- c()
  for(i in 1:number_layers) {
    new_cols <- c(new_cols, paste(cols, i, sep = "_"))
  }

  f <- new_cols

  print(data)
  if(nrow(data) <= number_layers) {
    row <- c()
    for(i in 1:nrow(data)) {
      if(na_omit & sum(is.na(data[i, ])) > 0) {
        next
      }
      row <- c(row, as.character(data[i, ]))
    }
    row <- c(row, rep(NA, length(f) - length(row)))

    f <- rbind(f, row)
  } else {
    for(i in 1:(length(data[, 1]) - number_layers + 1)) {
      time_slices <- data[i:(i + number_layers - 1), ]


      row <- c()
      for(j in 1:nrow(time_slices)) {

        row <- c(row, as.character(time_slices[j, ]))
      }

      f <- rbind(f, row)

    }
  }

  f <- as.data.frame(f)
  f <- f[-1, ]
  colnames(f) <- new_cols

  sorted_cols <- c()
  for(col in cols) {
    sorted_cols <- c(sorted_cols, paste(col, 1:number_layers, sep = "_"))
  }
  f <- f[, sorted_cols]
  return(f)
}


quantile_timeseries_data <- function(data, node_quantile) {
  data_quantile <- data
  cols <- colnames(data)
  breaks_list <- list()
  for(i in 1:ncol(data)) {

    values <- as.numeric(data[, i])
    value_max_min <- c(max(values), min(values))
    values <- scale_variables(values)

    breaks <- c(-Inf, unique(quantile(values, seq(0, 1, 1/node_quantile[i]))), Inf)
    breaks_list[[cols[i]]] <- list(breaks, value_max_min)
    data_quantile[, i] <- cut(values, breaks, labels = paste("Q", 1:(length(breaks) - 1), sep = ""))
  }

  for(i in 1:ncol(data)) {
    data_quantile[, i] <- factor(data_quantile[, i])
  }

  return(list(data_quantile, breaks_list))
}

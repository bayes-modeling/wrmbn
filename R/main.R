#' Preprocess time series data prepare for learning bayesian network
#'
#' @param data A data frame, each row is a time value and observations
#' @param type Bayesian network type, "discrete" or "continuous"
#' @param time_column Column name of "data", which values is time stamp
#' @param continuous_variable_names Column names of continuous variables
#' @param discrete_variable_names Column names of discrete variables
#' @param desire_layers Number layers of bayesian network, at least is 2
#' @param normalize_type Normalization type for continuous variables, "mean_normalization", "min_max" or "standardisation"
#' @param quantile_number Number of quantile level for type "discrete"
#' @param na_omit If true, NA/NaN values will be omit after preprocess
#' @return An object with pre-processsed data and preprocess parameters
#' @examples
#' library(wrmbn)
#'
#' data("data")
#' data("structure")
#'
#' uncorelate_node <- c("MTL", "QMX", "HCL", "CMB", "CTL", "CDX", "CBT")
#' for(node in uncorelate_node) {
#'  data <- data[, -which(colnames(data) == node)]
#' }
#'
#' for(node in uncorelate_node) {
#'   if(length(which(structure$from == node)) > 0) {
#'     structure <- structure[-which(structure$from == node), ]
#'   } else if(length(which(structure$to == node)) > 0) {
#'     structure <- structure[-which(structure$to == node), ]
#'   }
#' }

#' type <- "continuous"
#' time_column <- "date"
#' continuous_variable_names <- setdiff(colnames(data), time_column)
#' discrete_variable_names <- c()
#' desire_layers <- 3
#' normalize_type <- "mean_normalization"
#' preprocessed <- preprocess_timeseries_data(data, type, time_column,
#'                                            continuous_variable_names, discrete_variable_names, desire_layers,
#'                                            normalize_type, quantile_number = -1, na_omit = TRUE)
preprocess_timeseries_data <- function(data, type, time_column,
                                       continuous_variable_names, discrete_variable_names, desire_layers,
                                       normalize_type = NULL, quantile_number = -1, na_omit = TRUE) {

  # Read dataframe
  data <- data[, c(continuous_variable_names, discrete_variable_names)]

  # Reconstruct dataframe with desire layers
  data_constructed <- reconstruct_timeseries_data(data, desire_layers)
  data_constructed <- convert_data_type(data_constructed, continuous_variable_names, discrete_variable_names)
  continuous_variables <- c()
  for(name in continuous_variable_names) {
    continuous_variables <- c(continuous_variables, paste(name, 1:desire_layers, sep = "_"))
  }

  discrete_variables <- c()
  for(name in discrete_variable_names) {
    discrete_variables <- c(discrete_variables, paste(name, 1:desire_layers, sep = "_"))
  }

  data_constructed <- data_constructed[, c(continuous_variables, discrete_variables)]
  #Normalize data
  normalizers <- NULL
  if(!is.null(normalize_type)) {
    data_normalized <- normalize_data(data_constructed, normalize_type)
    data_constructed <- data_normalized$data
    normalizers <- data_normalized$normalizers
  }

  result <- list()
  result[["continuous_variables"]] <- continuous_variables
  result[["discrete_variables"]] <- discrete_variables
  result[["quantile_number"]] <- quantile_number
  result[["desire_layers"]] <- desire_layers
  result[["time_column"]] <- time_column
  result[["type"]] <- type
  result[["normalizers"]] <- normalizers
  result[["normalize_tye"]] <- normalize_type
  if(type == "discrete") {
    quantiled_variables <- c()

    if(quantile_number > 1) {
      quantile_number <- rep(quantile_number, length(continuous_variable_names) * desire_layers)
      for(variable in continuous_variable_names) {
        quantiled_variables <- c(quantiled_variables, paste(variable, 1:desire_layers, sep = "_"))
      }
    }

    data_discrete_dbn <- quantile_timeseries_data(data_constructed[, quantiled_variables], quantile_number)
    break_list <- data_discrete_dbn[[2]]
    data_discrete <- data_discrete_dbn[[1]]

    result[["data"]] <- data_discrete
    result[["break_list"]] <- break_list
  } else {
    result[["data"]] <- data_constructed
  }

  return(result)
}


#' Preprocess time series data prepare for test learned bayesian network
#'
#' @param data A data frame, each row is a time value and observations
#' @param continuous_variable_names Column names of continuous variables
#' @param discrete_variable_names Column names of discrete variables
#' @param desire_layers Number layers of bayesian network, at least is 2
#' @param time_column Column name of "data", which values is time stamp, default is NULL
#' @param normalize_type Normalization type for continuous variables, "mean_normalization", "min_max" or "standardisation", default is NULL
#' @param normalizers Normalize parameters, default is NULL
#' @return An list object, each element is processed data frame
#' @examples
#' library(wrmbn)
#' data("preprocessed")
#' head(data)
#' continuous_variables <- preprocessed$continuous_variables
#' discrete_variables <- preprocessed$discrete_variables
#' desire_layers <- preprocessed$desire_layers
#' time_column <- "date"
#' normalize_type <- preprocessed$normalize_tye
#' normalizers <- preprocessed$normalizers
#' test_data <- preprocess_test_data(data, continuous_variables, discrete_variables, desire_layers,
#'                                   time_column, normalize_type, normalizers)
preprocess_test_data <- function(data, continuous_variables, discrete_variables, desire_layers,
                                 time_column = NULL, normalize_type = NULL, normalizers = NULL) {
  time_values <- NULL
  if(!is.null(time_column)) {
    time_values <- data[, time_column]
  }
  data_constructed <- reconstruct_timeseries_data(data, desire_layers, FALSE)
  data_constructed <- data_constructed[, c(continuous_variables, discrete_variables)]
  data_constructed <- convert_data_type(data_constructed, continuous_variables, discrete_variables)

  if(!is.null(normalize_type) & !is.null(normalizers)) {
    data_constructed <- normalize_data(data_constructed, normalize_type, normalizers)$data
  }

  return(list(time_values=time_values, data=data_constructed))
}

#' Filter blacklist and whitelist before training bayesian network
#'
#' @param continuous_data A data frame which rows are financial reports and columns is KPI pre-processed
#' @param desire_layers Number of layers for training bayesian network, less than or equal current layers
#' @param quantile_number For type is "discrete", number of quantile for each continuous variable
#' @param continuous_variables Column names of continuous dynamic KPI
#' @param discrete_variables Column names of discrete dynamic KPI
#' @param known_structure Relation between KPI specific by user
#' @param corr_threshold Threshold for filter arcs with low correlation
#' @param is_blacklist_internal Is blacklist all internal arcs in n-i layers (n is number layers, i = 1..n-1)
#' @param is_variable_only Is known_structure is KPI only, without time notation
#' @param is_blacklist_other Is blacklist all except whitelist
#' @param custom_blacklist Blacklist specific by user
#' @param custom_whitelist Whitelist sepcific by user
#' @return An object with blacklist and whitelist of network
#' @examples
#' data("preprocessed")
#' data("data_small")
#' continuous_data <- preprocessed$data
#' desire_layers <- preprocessed$desire_layers
#' quantile_number <- preprocessed$quantile_number
#' continuous_dynamic_variables <- preprocessed$continuous_dynamic_variables
#' continuous_static_variables <- preprocessed$continuous_static_variables
#' discrete_static_variables <- preprocessed$discrete_static_variables
#' known_structure <- NULL
#' corr_threshold <- 0.4
#' is_blacklist_internal <- TRUE
#' is_variable_only <- TRUE
#' is_blacklist_other <- FALSE
#' custom_blacklist <- NULL
#' custom_whitelist <- NULL

#' bl_wl <- get_continuous_structure_filter(continuous_data, desire_layers, quantile_number, continuous_dynamic_variables,
#'                                                 continuous_static_variables, discrete_static_variables,
#'                                                 known_structure, corr_threshold, is_blacklist_internal,
#'                                                 is_variable_only, is_blacklist_other,
#'                                                 custom_blacklist, custom_whitelist)
get_continuous_structure_filter <- function(continuous_data, desire_layers, quantile_number,
                                            continuous_variables, discrete_variables,
                                            known_structure, corr_threshold, is_blacklist_internal,
                                            is_variable_only, is_blacklist_other,
                                            custom_blacklist, custom_whitelist){


  continuous_data <- continuous_data[, c(continuous_variables, discrete_variables)]

  rho <- abs(cor(convert_variables_to_factor(continuous_data)$data))
  rho <- as.matrix(rho)
  arc_structures <- structure_filter(known_structure, desire_layers,
                                     c(continuous_variables, discrete_variables),
                                     is_blacklist_internal, is_blacklist_other, is_variable_only,
                                     rho, corr_threshold,
                                     custom_blacklist, custom_whitelist)

  return(arc_structures)
}

#' Filter blacklist and whitelist before training bayesian network
#'
#' @param continuous_data A data frame which rows are financial reports and columns is KPI pre-processed
#' @param desire_layers Number of layers for training bayesian network, less than or equal current layers
#' @param quantile_number For type is "discrete", number of quantile for each continuous variable
#' @param continuous_variables Column names of continuous dynamic KPI
#' @param discrete_variables Column names of discrete dynamic KPI
#' @param known_structure Relation between KPI specific by user
#' @param corr_threshold Threshold for filter arcs with low correlation
#' @param is_blacklist_internal Is blacklist all internal arcs in n-i layers (n is number layers, i = 1..n-1)
#' @param is_variable_only Is known_structure is KPI only, without time notation
#' @param is_blacklist_other Is blacklist all except whitelist
#' @param custom_blacklist Blacklist specific by user
#' @param custom_whitelist Whitelist sepcific by user
#' @return An object with blacklist and whitelist of network
#' @examples
#' data("preprocessed")
#' data("data_small")
#' continuous_data <- preprocessed$data
#' desire_layers <- preprocessed$desire_layers
#' quantile_number <- preprocessed$quantile_number
#' continuous_dynamic_variables <- preprocessed$continuous_dynamic_variables
#' continuous_static_variables <- preprocessed$continuous_static_variables
#' discrete_static_variables <- preprocessed$discrete_static_variables
#' known_structure <- NULL
#' corr_threshold <- 0.4
#' is_blacklist_internal <- TRUE
#' is_variable_only <- TRUE
#' is_blacklist_other <- FALSE
#' custom_blacklist <- NULL
#' custom_whitelist <- NULL

#' bl_wl <- tfdbn::get_continuous_structure_filter(continuous_data, desire_layers, quantile_number, continuous_dynamic_variables,
#'                                                 continuous_static_variables, discrete_static_variables,
#'                                                 known_structure, corr_threshold, is_blacklist_internal,
#'                                                 is_variable_only, is_blacklist_other,
#'                                                 custom_blacklist, custom_whitelist)
get_discrete_structure_filter <- function(discrete_data, desire_layers, quantile_number,
                                          continuous_variables, discrete_variables,
                                          known_structure, corr_threshold, is_blacklist_internal,
                                          is_variable_only, is_blacklist_other,
                                          custom_blacklist, custom_whitelist, debug = FALSE){


  discrete_data <- discrete_data[, c(continuous_variables, discrete_variables)]

  rho <- abs(cor(convert_variables_to_factor(discrete_data)$data))
  rho <- as.matrix(rho)

  arc_structures <- structure_filter(known_structure, desire_layers,
                                     c(continuous_variables, discrete_variables),
                                     is_blacklist_internal, is_blacklist_other, is_variable_only,
                                     rho, corr_threshold,
                                     custom_blacklist, custom_whitelist)

  return(arc_structures)
}

#' Training bayesian network
#'
#' @param training_type Training type, discrete or continuous
#' @param data Data for training
#' @param number_layers Number of bayesian network layer
#' @param bl List of blacklist
#' @param wl List of whitelist
#' @param n_cluster Number of core for training
#' @param algorithms Algorithms for learning network structure and parameters
#' @param number_bootstrap Number of bootstrap
#' @param debug Debug mode
#' @return An object list of trained model and training parameters
#' @examples
#' library(wrmbn)
#' data("preprocessed")
#' training_type <- preprocessed$type
#' data <- preprocessed$data
#' number_layers <- preprocessed$desire_layers
#' bl <- bl_wl$blacklist
#' wl <- bl_wl$whitelist
#' n_cluster <- 4
#' algorithms <- c("gs", "hc", "tabu", "iamb", "inter.iamb", "fast.iamb")
#' algorithms <- c("hc", "tabu")
#' trained_models <- training_model(training_type, data, number_layers, bl, wl, n_cluster, algorithms, number_bootstrap = 500, debug = FALSE)
training_model <- function(training_type, data, number_layers, bl, wl, n_cluster, algorithms, number_bootstrap = 100, debug = FALSE) {
  if(debug) {
    cat("Training model ", training_type, "\n")
  }

  cluster <- parallel::makeCluster(n_cluster)

  data_converted <- convert_variables_to_factor(data)

  data_factors <- data_converted$data_factors
  data <- data_converted$data

  bl <- remove_unknown_arc(bl, colnames(data))
  wl <- remove_unknown_arc(wl, colnames(data))
  results <- list()
  results$data_factors <- data_factors
  for(algor in algorithms) {

    if(debug) {
      cat("Training model by algorithm ", algor, "\n")
    }

    begin <- Sys.time()
    trained <- structure_learning(data, desire_layers, bl, wl, slearning_algo = algor, number_bootstrap = number_bootstrap, cluster, debug)
    end <- Sys.time()
    trained$time <- begin - end
    trained$date <- date()
    trained$blacklist <- bl
    trained$whitelist <- wl

    if(debug) {
      cat("Training model by algorithm ", algor, "compeleted", "\n")
    }

    results[[algor]] <- trained

  }

  parallel::stopCluster(cluster)

  if(debug) {
    cat("Training model", training_type, "compeleted", "\n")
  }


  return(results)
}

#' Cross validation for learning bayesian network
#'
#' @param data Processed data frame for learning bayesian network
#' @param algorithms Algorithms for cross validation
#' @param wl Network white list
#' @param bl Network black list
#' @param target A character string, the label of target node for prediction
#' @param k_fold The data are split in k subsets of equal size
#' @param runs A positive integer number, the number of times k-fold or hold-out cross-validation will be run.
#' @param loss A character string, the label of a loss function, detail see https://cran.r-project.org/web/packages/bnlearn/bnlearn.pdf
#' @param n_cluster an optional cluster object from package parallel.
#' @param debug Debug mode
#' @return Cross validation result for each algorithm
#' @example
#' library(wrmbn)
#' data("preprocessed")
#' data <- preprocessed$data
#' algorithms <- c("gs", "hc", "tabu", "iamb", "inter.iamb", "fast.iamb")
#' algorithms <- c("hc", "tabu")
#' target <- "MBT_3"
#' wl <- bl_wl$whitelist
#' bl <- bl_wl$blacklist

#' result <- wrmbn::cross_validation_learning_algorithms(data, algorithms, wl, bl, target, n_cluster = 4, debug = TRUE)
#' boxplot(result$hc, result$tabu, names = algorithms)
cross_validation_learning_algorithms <- function(data, algorithms, wl, bl, target,
                                                 k_fold = 5, runs = 10, loss = "mse",
                                                 n_cluster = NULL, debug = NULL) {
  if(debug) {
    cat("Cross validation for algorithms", paste(algorithms, collapse = ","), "with k_fold =", k_fold, "and target =", target, "\n")
  }
  cluster <- NULL
  if(is.null(n_cluster)) {
    cluster <- parallel::makeCluster(n_cluster)
  }

  result <- list()

  for(algorithm in algorithms) {
    if(debug) {
      cat("Cross validation for algorithm", algorithm, "with k_fold =", k_fold, "and target =", target, "\n")
    }

    algorithm.cv <- bnlearn::bn.cv(data, k = k_fold, algorithm, loss = loss, runs = runs,
                                   algorithm.args = list(whitelist = wl, blacklist = bl),
                                   loss.args = list(target = target),
                                   cluster = cluster)

    result[[algorithm]] <- extract_cross_validation_loss(algorithm.cv)

    if(debug) {
      cat("Cross validation for algorithm", algorithm, "with k_fold =", k_fold, "and target =", target, "completed", "\n")
    }
  }

  if(!is.null(cluster)) {
    parallel::stopCluster(cluster)
  }

  if(debug) {
    cat("Cross validation for algorithms", paste(algorithms, collapse = ","), "with k_fold =", k_fold, "and target =", target, "completed", "\n")
  }

  return(result)
}

#' Impute missing data in data frame
#'
#' @param fitted A fitted model
#' @param data A data frame, contains missing data, each row is a time value and observations, structure is the same as learning data before process
#' @param continuous_variables Column names of continuous variables
#' @param discrete_variables Column names of discrete variables
#' @param time_column Column name of "data", which values is time stamp, default is NULL
#' @param normalize_type Normalization type for continuous variables, "mean_normalization", "min_max" or "standardisation", default is NULL
#' @param normalizers Normalize parameters, default is NULL
#' @return A full filled data frame
#' @examples
#' library(wrmbn)
#' data("data")
#' data("preprocessed")
#' data("trained_models")
#' head(data)
#' continuous_variables <- preprocessed$continuous_variables
#' discrete_variables <- preprocessed$discrete_variables
#' desire_layers <- preprocessed$desire_layers
#' time_column <- "date"
#' normalize_type <- preprocessed$normalize_tye
#' normalizers <- preprocessed$normalizers
#' fitted <- trained_models$hc$fitted
#'
#' imputed_data <- impute_missing_data(fitted, data, continuous_variables, discrete_variables,
#'                                     time_column, normalize_type, normalizers, debug = FALSE)
#' head(imputed_data)
impute_missing_data <- function(fitted, data, continuous_variables, discrete_variables,
                                time_column = NULL, normalize_type = NULL, normalizers = NULL, debug = FALSE) {
  if(debug) {
    cat("Impute missing for data with dim", dim(data))
  }

  #Reconstruct and normalize data (if need)
  test_data <- preprocess_test_data(data, continuous_variables, discrete_variables, desire_layers,
                                           time_column, normalize_type, normalizers)

  time_values <- test_data$time_values
  processed_data <- test_data$data

  #Impute missing data
  impute_data <- bnlearn::impute(fitted, processed_data)

  #Re-scale data to actual
  impute_data <- reverse_normalized_data(impute_data, normalize_type, normalizers)

  #Reverse constructed data
  impute_data <- reverse_constructed_data(impute_data, desire_layers)
  all_names <- colnames(impute_data)

  if(!is.null(time_column)) {
    impute_data <- impute_data[1:length(time_values), ]
    impute_data <- cbind(time_values, impute_data)
    colnames(impute_data) <- c(time_column, all_names)
  }


  if(debug) {
    cat("Impute missing for data with dim", dim(data), "completed")
  }

  return(impute_data)
}

#' Predict data
#'
#' @param fitted A fitted model
#' @param predictor_data A data frame, each row is evidence
#' @param predicted_nodes List of variables need to predict
#' @param number_layers Number layers of fitted model
#' @param continuous_variables Column names of continuous variables
#' @param discrete_variables Column names of discrete variables
#' @param time_column Column name of "data", which values is time stamp, default is NULL
#' @param time_format Time format of time column, default is %m/%d/%y
#' @param normalize_type Normalization type for continuous variables, "mean_normalization", "min_max" or "standardisation", default is NULL
#' @param normalizers Normalize parameters, default is NULL
#' @return A data frame with evidence and predicted variables
#' @examples
#' library(wrmbn)
#' data("data")
#' data("preprocessed")
#' data("trained_models")
#' predictor_data <- data[2019:2020, c("date", "HND", "HCT")]
#' print(predictor_data)
#' normalizers <- preprocessed$normalizers
#' normalize_type <- preprocessed$normalize_tye
#' fitted <- trained_models$hc$fitted
#' continuous_variables <- c("HND", "HCT")
#' discrete_variables <- c()
#' predicted_nodes <- c("MBT")
#' number_layers <- 3
#'
#' predicted_data <- predict_data(fitted, predictor_data, predicted_nodes, number_layers,
#'                                continuous_variables, discrete_variables,
#'                                time_column = "date", "%m/%d/%y", normalizers, normalize_type,
#'                                method = "lw", TRUE)
#'
#' print(predicted_data)
#' actual_data <- data[2019:2021, c("date", "HND", "HCT", "MBT")]
#' print(actual_data)
predict_data <- function(fitted, predictor_data, predicted_nodes, number_layers,
                         continuous_variables, discrete_variables,
                         time_column = NULL, time_format = "%m/%d/%y", normalizers = NULL, normalize_type = NULL,
                         method = "lw", debug = FALSE) {
  if(debug) {
    cat("Predict for nodes", predicted_nodes, "\n")
  }
  predictor_data <- predictor_data[, c(time_column, continuous_variables, discrete_variables)]

  #Bind predicted values
  for(node in predicted_nodes) {
    predictor_data <- cbind(predictor_data, rep(NA, nrow(predictor_data)))
  }

  colnames(predictor_data) <- c(time_column, continuous_variables, discrete_variables, predicted_nodes)

  predictor_data_without_time <- predictor_data
  time_values <- NULL
  if(!is.null(time_column)) {
    time_values <- as.Date(predictor_data[, time_column], time_format)
    predictor_data_without_time <- predictor_data[, -which(colnames(predictor_data) == time_column)]
  }

  continuous_variable_names <- generate_variables(continuous_variables, number_layers)
  data_constructed <- reconstruct_timeseries_data(predictor_data_without_time, number_layers, FALSE)
  data_constructed <- convert_data_type(data_constructed, c(continuous_variables, predicted_nodes), c())
  if(!is.null(normalize_type) & !is.null(normalizers)) {
    data_constructed <- normalize_data(data_constructed, normalize_type, normalizers)$data
  }

  #Predict values
  for(i in seq(nrow(data_constructed))) {
    while(TRUE) {
      predictor <- data_constructed[i, ]
      na_columns <- c()
      for(col in colnames(predictor)) {
        if(sum(is.na(predictor[, col])) > 0) {
          na_columns <- c(na_columns, col)
        }
      }
      evidence_columns <- setdiff(colnames(predictor), na_columns)
      if(length(na_columns) == 0) {
        break
      }

      predict_node <- na_columns[1]
      predictor_node <- data_constructed[, evidence_columns]
      if(debug) {
        cat("Predict for node", predict_node, "with predictors", evidence_columns, "\n")
      }
      dist <- bnlearn::cpdist(fitted, nodes = predict_node,
                              evidence = as.list(predictor_node), method = method)
      predicted_value = mean(dist[, predict_node])
      data_constructed[i, predict_node] <- predicted_value

      if(debug) {
        cat("Predict for node", predict_node, "with predictors", evidence_columns, "completed with value", predicted_value, "\n")
      }
    }

  }

  data_constructed <- reverse_normalized_data(data_constructed, normalize_type, normalizers)
  data_constructed <- reverse_constructed_data(data_constructed, number_layers)

  if(debug) {
    cat("Predict for nodes", predicted_nodes, "completed", "\n")
  }

  if(!is.null(time_column) & !is.null(time_values)) {

    actual_length = length(time_values)
    predicted_length = nrow(data_constructed)
    new_time_values <- NULL
    if(actual_length < predicted_length) {
      new_time_values <- seq.Date(time_values[1], length = predicted_length, by = "day")
      data_constructed <- cbind(new_time_values, data_constructed)
    } else {
      data_constructed <- cbind(time_values, data_constructed)
    }

    colnames(data_constructed) <- c(time_column, continuous_variables, discrete_variables, predicted_nodes)
  }
  return(data_constructed)
}





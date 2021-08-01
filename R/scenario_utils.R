
dist_scenario <- function(fitted, data, event_nodes, evidence_nodes, type, range, method = "lw") {
  evidence_equation <- ""
  nodes <- generate_dist_node(event_nodes)
  result <- list()
  for(i in 1:nrow(data)) {
    if(type == "discrete") {
      evidence_equation <- generate_discrete_equation(data[1, ], evidence_nodes, "evidence")
      dist <- paste(paste("cpdist(fitted, ", nodes, ", ", evidence_equation, ")", sep = ""))
      result[[i]] <- eval(parse(text=dist))
    } else {
      predictor <- data[i, evidence_nodes]
      print(event_nodes)
      dist <- cpdist(fitted, nodes = event_nodes,
                     evidence = as.list(predictor), method = method)
      result[[i]] <- dist
    }
  }
  
  return(result)
}

query_scenario <- function(fitted, data, event_nodes, evidence_nodes, type, range) {
  query <- ""
  result <- c()
  for(i in 1:nrow(data)) {
    if(type == "discrete") {
      
      evidence_equation <- generate_discrete_equation(data[i, ], evidence_nodes, "evidence")
      event_equation <- generate_discrete_equation(data[i, ], event_nodes, "event")
      query <- paste(paste("cpquery(fitted, ", event_equation, ", ", evidence_equation, ", n = 1e6)", sep = ""))
      
    } else {
      
      evidence_equation <- generate_continuous_equation(data[i, ], evidence_nodes, "evidence", range)
      event_equation <- generate_continuous_equation(data[i, ], event_nodes, "event", range)
      query <- paste(paste("cpquery(fitted, ", event_equation, ", ", evidence_equation, ", n = 1e6)", sep = ""))
    }
    
    print(query)
    
    result <- c(result, eval(parse(text=query)))
  }
    
  return(result)
}

query_sector_kpi <- function(fitted, query_variables, value, type) {
  
  sector_dist <- cpdist(fitted, nodes = query_variables,
                 evidence = (NGANHKT == value), method = method)
  
  other_dist <- cpdist(fitted, nodes = query_variables,
                        evidence = (NGANHKT != value), method = method)
  
  
}

generate_discrete_equation <- function(data, nodes, equation_type) {
  equation <- c()
  for(i in 1:length(nodes)) {
    event <- nodes[i]
    for(j in 1:nrow(data)) {
      value <- data[j, event]
      event_string <- as.character(event)
      equation = c(equation, paste(event_string, " == \"", value, "\"", sep = ""))
    }
  }
  
  equation <- paste(equation_type, " = (", paste(equation, collapse = " & "), ")", sep = "")
  
  return(equation)
}

generate_continuous_equation <- function(data, nodes, equation_type, range) {
  equation <- c()
  for(i in 1:length(nodes)) {
    for(j in 1:nrow(data)) {
      event <- nodes[i]
      value <- data[j, event]
      value_upper <- value + range
      value_lower <- value - range
      
      event_string <- as.character(event)
      
      equation = c(equation, paste("(", event_string, " >= ", value_lower, " & ", event_string, " <= ", value_upper, ")", sep = ""))
    }
    
    
  }
  
  print(equation)
  
  equation <- paste(equation_type, " = (", paste(equation, collapse = " & "), ")", sep = "")
  
  print(equation)
  
  return(equation)
}

generate_dist_node <- function(nodes) {
  node_strings <- c()
  for(node in nodes) {
    node_strings <- c(node_strings, paste("\"", node, "\"", sep = ""))
  }
  equation <- paste("nodes = c(", paste(node_strings, collapse = ","), ")", sep = "")
  print(equation)
  return(equation)
}


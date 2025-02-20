#' Sort population by fitness
#'
#' @param population List of individuals with fitness values
#' @param minimize Whether to minimize (TRUE) or maximize (FALSE) fitness
#' @return Sorted population
#' @keywords internal
sort_by_fitness <- function(population, minimize = TRUE) {
  if (minimize) {
    population[order(sapply(population, function(ind) ind$fitness))]
  } else {
    population[order(sapply(population, function(ind) ind$fitness), decreasing = TRUE)]
  }
}

#' Validate that a genotype conforms to grammar structure
#'
#' @param genotype Genotype to validate
#' @param grammar Grammar structure
#' @return TRUE if valid, FALSE otherwise
#' @export
validate_genotype <- function(genotype, grammar) {
  # Check that genotype is a list
  if (!is.list(genotype)) {
    return(FALSE)
  }
  
  # All values should be numeric between 0 and 1
  for (nt in names(genotype)) {
    if (!(nt %in% grammar$non_terminals)) {
      return(FALSE)
    }
    
    if (!is.numeric(genotype[[nt]]) || 
        any(genotype[[nt]] < 0) || 
        any(genotype[[nt]] > 1)) {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Visualize the evolution history
#'
#' @param history Evolution history data frame
#' @param title Plot title
#' @return ggplot object
#' @export
plot_evolution <- function(history, title = "PSGE Evolution") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  p <- ggplot2::ggplot(history, ggplot2::aes(x = generation)) +
    ggplot2::geom_line(ggplot2::aes(y = best_fitness, color = "Best Fitness")) +
    ggplot2::geom_line(ggplot2::aes(y = mean_fitness, color = "Mean Fitness")) +
    ggplot2::scale_color_manual(values = c("Best Fitness" = "blue", "Mean Fitness" = "red")) +
    ggplot2::labs(
      title = title,
      x = "Generation",
      y = "Fitness",
      color = "Metric"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}

#' Convert individual to string representation
#'
#' @param individual Individual to convert
#' @param include_genotype Whether to include genotype details
#' @return String representation of individual
#' @export
individual_to_string <- function(individual, include_genotype = FALSE) {
  result <- paste0(
    "Fitness: ", round(individual$fitness, 6), "\n",
    "Phenotype: ", individual$phenotype, "\n",
    "Tree depth: ", individual$tree_depth, "\n"
  )
  
  if (include_genotype) {
    result <- paste0(
      result,
      "Genotype:\n"
    )
    
    for (nt in names(individual$genotype)) {
      result <- paste0(
        result,
        "  ", nt, ": [", 
        paste(round(individual$genotype[[nt]], 4), collapse = ", "),
        "]\n"
      )
    }
  }
  
  return(result)
}

#' Serialize grammar to JSON
#'
#' @param grammar Grammar structure to serialize
#' @param file Output file path (optional)
#' @return JSON string if file is NULL, otherwise writes to file
#' @export
grammar_to_json <- function(grammar, file = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  json <- jsonlite::toJSON(grammar, pretty = TRUE, auto_unbox = TRUE)
  
  if (!is.null(file)) {
    writeLines(json, file)
    invisible(json)
  } else {
    return(json)
  }
}

#' Load grammar from JSON
#'
#' @param json_str JSON string or file path
#' @param is_file Whether json_str is a file path
#' @return Grammar structure
#' @export
grammar_from_json <- function(json_str, is_file = FALSE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (is_file) {
    if (!file.exists(json_str)) {
      stop("File not found: ", json_str)
    }
    json_str <- paste(readLines(json_str), collapse = "\n")
  }
  
  grammar <- jsonlite::fromJSON(json_str, simplifyVector = FALSE)
  
  # Validate the loaded grammar
  grammar <- validate_grammar(grammar)
  
  return(grammar)
}

#' Generate a representative sample of phenotypes from a grammar
#'
#' @param grammar Grammar structure
#' @param n Number of samples to generate
#' @param max_depth Maximum tree depth
#' @return List of generated phenotypes
#' @export
sample_phenotypes <- function(grammar, n = 10, max_depth = 10) {
  phenotypes <- character(n)
  
  for (i in seq_len(n)) {
    # Generate random individual
    genotype <- generate_random_individual(grammar, max_depth)
    
    # Map to phenotype
    mapping_result <- map_genotype(grammar, genotype, max_depth)
    phenotypes[i] <- mapping_result$phenotype
  }
  
  return(phenotypes)
}

#' Export a grammar back to BNF format
#'
#' @param grammar Grammar structure
#' @param file Output file path (optional)
#' @param include_probs Whether to include probabilities as comments
#' @return BNF string if file is NULL, otherwise writes to file
#' @export
grammar_to_bnf <- function(grammar, file = NULL, include_probs = TRUE) {
  lines <- character(0)
  
  for (nt in grammar$non_terminals) {
    line <- paste0("<", nt, "> ::= ")
    rules <- grammar$rules[[nt]]
    
    for (i in seq_along(rules)) {
      if (i > 1) {
        line <- paste0(line, " | ")
      }
      
      # Format symbols
      symbols <- character(length(rules[[i]]$symbols))
      for (j in seq_along(rules[[i]]$symbols)) {
        if (rules[[i]]$types[j] == "NT") {
          symbols[j] <- paste0("<", rules[[i]]$symbols[j], ">")
        } else {
          symbols[j] <- rules[[i]]$symbols[j]
        }
      }
      
      line <- paste0(line, paste(symbols, collapse = " "))
      
      # Add probability as comment if requested
      if (include_probs) {
        line <- paste0(line, " # [", format(rules[[i]]$prob, digits = 4), "]")
      }
    }
    
    lines <- c(lines, line)
  }
  
  bnf <- paste(lines, collapse = "\n")
  
  if (!is.null(file)) {
    writeLines(bnf, file)
    invisible(bnf)
  } else {
    return(bnf)
  }
}

#' Calculate statistics about grammar usage
#'
#' @param grammar Grammar structure
#' @param population Population of individuals
#' @return Data frame with statistics
#' @export
analyze_grammar_usage <- function(grammar, population) {
  # Initialize counters
  rule_usage <- list()
  for (nt in grammar$non_terminals) {
    rule_usage[[nt]] <- rep(0, length(grammar$rules[[nt]]))
  }
  
  # Count rule usage across population
  for (ind in population) {
    if (is.null(ind$genotype)) {
      next
    }
    
    for (nt in names(ind$genotype)) {
      if (!(nt %in% names(grammar$rules))) {
        next
      }
      
      for (codon in ind$genotype[[nt]]) {
        # Determine which rule this codon selected
        cum_prob <- 0
        for (i in seq_along(grammar$rules[[nt]])) {
          cum_prob <- cum_prob + grammar$rules[[nt]][[i]]$prob
          if (codon <= cum_prob) {
            rule_usage[[nt]][i] <- rule_usage[[nt]][i] + 1
            break
          }
        }
      }
    }
  }
  
  # Convert to data frame
  result <- data.frame(
    non_terminal = character(0),
    rule_index = integer(0),
    rule_repr = character(0),
    usage_count = integer(0),
    usage_percent = numeric(0),
    current_prob = numeric(0)
  )
  
  for (nt in names(rule_usage)) {
    total_usage <- sum(rule_usage[[nt]])
    
    for (i in seq_along(rule_usage[[nt]])) {
      # Create rule representation
      rule_symbols <- character(length(grammar$rules[[nt]][[i]]$symbols))
      for (j in seq_along(grammar$rules[[nt]][[i]]$symbols)) {
        if (grammar$rules[[nt]][[i]]$types[j] == "NT") {
          rule_symbols[j] <- paste0("<", grammar$rules[[nt]][[i]]$symbols[j], ">")
        } else {
          rule_symbols[j] <- grammar$rules[[nt]][[i]]$symbols[j]
        }
      }
      rule_repr <- paste(rule_symbols, collapse = " ")
      
      # Calculate percentage
      usage_percent <- if (total_usage > 0) {
        100 * rule_usage[[nt]][i] / total_usage
      } else {
        0
      }
      
      # Add to result
      result <- rbind(result, data.frame(
        non_terminal = nt,
        rule_index = i,
        rule_repr = rule_repr,
        usage_count = rule_usage[[nt]][i],
        usage_percent = usage_percent,
        current_prob = 100 * grammar$rules[[nt]][[i]]$prob
      ))
    }
  }
  
  return(result)
}

#' Find minimum tree depth required for each non-terminal
#'
#' @param grammar Grammar structure
#' @return Named vector with minimum depths
#' @export
calculate_min_depths <- function(grammar) {
  # Initialize depths to Inf
  depths <- setNames(rep(Inf, length(grammar$non_terminals)), grammar$non_terminals)
  
  # Terminals have depth 0
  terminal_depths <- setNames(rep(0, length(grammar$terminals)), grammar$terminals)
  
  # Keep iterating until no changes
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    
    for (nt in grammar$non_terminals) {
      for (rule in grammar$rules[[nt]]) {
        # Calculate depth for this rule
        max_symbol_depth <- 0
        all_symbols_known <- TRUE
        
        for (i in seq_along(rule$symbols)) {
          if (rule$types[i] == "NT") {
            if (is.infinite(depths[rule$symbols[i]])) {
              all_symbols_known <- FALSE
              break
            }
            max_symbol_depth <- max(max_symbol_depth, depths[rule$symbols[i]])
          }
          # Terminals have depth 0, so they don't increase max_symbol_depth
        }
        
        if (all_symbols_known) {
          rule_depth <- max_symbol_depth + 1
          if (rule_depth < depths[nt]) {
            depths[nt] <- rule_depth
            changed <- TRUE
          }
        }
      }
    }
  }
  
  return(depths)
}

#' Create a visualization of grammar derivations
#'
#' @param grammar Grammar structure
#' @param start_symbol Starting non-terminal (defaults to grammar's start symbol)
#' @param max_depth Maximum depth for visualization
#' @param show_probs Whether to show probabilities
#' @return DiagrammeR graph object
#' @export
visualize_grammar <- function(grammar, start_symbol = NULL, max_depth = 3, show_probs = TRUE) {
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) {
    stop("Package 'DiagrammeR' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (is.null(start_symbol)) {
    start_symbol <- grammar$start_symbol
  }
  
  # Initialize nodes and edges
  nodes <- data.frame(
    id = 1,
    label = paste0("<", start_symbol, ">"),
    type = "NT",
    stringsAsFactors = FALSE
  )
  
  edges <- data.frame(
    from = integer(0),
    to = integer(0),
    label = character(0),
    stringsAsFactors = FALSE
  )
  
  # Counter for node IDs
  next_id <- 2
  
  # Queue for BFS
  queue <- list(list(symbol = start_symbol, node_id = 1, depth = 1))
  
  # BFS traversal
  while (length(queue) > 0) {
    current <- queue[[1]]
    queue <- queue[-1]
    
    if (current$depth > max_depth) {
      next
    }
    
    # Process rules for this non-terminal
    if (current$symbol %in% names(grammar$rules)) {
      for (i in seq_along(grammar$rules[[current$symbol]])) {
        rule <- grammar$rules[[current$symbol]][[i]]
        
        # Create label for the rule
        rule_label <- if (show_probs) {
          paste0("Rule ", i, " [", round(rule$prob, 2), "]")
        } else {
          paste0("Rule ", i)
        }
        
        # Add node for this rule
        rule_node_id <- next_id
        next_id <- next_id + 1
        
        nodes <- rbind(nodes, data.frame(
          id = rule_node_id,
          label = rule_label,
          type = "rule",
          stringsAsFactors = FALSE
        ))
        
        # Connect parent node to rule node
        edges <- rbind(edges, data.frame(
          from = current$node_id,
          to = rule_node_id,
          label = "",
          stringsAsFactors = FALSE
        ))
        
        # Process symbols in the rule
        for (j in seq_along(rule$symbols)) {
          symbol_type <- rule$types[j]
          symbol_value <- rule$symbols[j]
          
          # Format symbol label
          symbol_label <- if (symbol_type == "NT") {
            paste0("<", symbol_value, ">")
          } else {
            paste0('"', symbol_value, '"')
          }
          
          # Add node for this symbol
          symbol_node_id <- next_id
          next_id <- next_id + 1
          
          nodes <- rbind(nodes, data.frame(
            id = symbol_node_id,
            label = symbol_label,
            type = symbol_type,
            stringsAsFactors = FALSE
          ))
          
          # Connect rule node to symbol node
          edges <- rbind(edges, data.frame(
            from = rule_node_id,
            to = symbol_node_id,
            label = paste0("pos ", j),
            stringsAsFactors = FALSE
          ))
          
          # Add non-terminal to queue for further expansion
          if (symbol_type == "NT" && current$depth < max_depth) {
            queue <- c(queue, list(list(
              symbol = symbol_value,
              node_id = symbol_node_id,
              depth = current$depth + 1
            )))
          }
        }
      }
    }
  }
  
  # Create graph
  graph <- DiagrammeR::create_graph() %>%
    DiagrammeR::add_nodes_from_table(
      table = nodes,
      label_col = "label",
      type_col = "type"
    ) %>%
    DiagrammeR::add_edges_from_table(
      table = edges,
      from_col = "from",
      to_col = "to",
      rel_col = "label"
    ) %>%
    DiagrammeR::add_global_graph_attrs(
      attr_type = "graph",
      attr = "layout",
      value = "dot"
    ) %>%
    DiagrammeR::add_global_graph_attrs(
      attr_type = "graph",
      attr = "rankdir",
      value = "TB"
    )
  
  # Style nodes by type
  graph <- graph %>%
    DiagrammeR::node_aes(
      node_attr = "style",
      values = "filled",
      node_sel = nodes$type == "NT"
    ) %>%
    DiagrammeR::node_aes(
      node_attr = "fillcolor",
      values = "lightblue",
      node_sel = nodes$type == "NT"
    ) %>%
    DiagrammeR::node_aes(
      node_attr = "shape",
      values = "box",
      node_sel = nodes$type == "rule"
    ) %>%
    DiagrammeR::node_aes(
      node_attr = "style",
      values = "filled",
      node_sel = nodes$type == "T"
    ) %>%
    DiagrammeR::node_aes(
      node_attr = "fillcolor",
      values = "lightgreen",
      node_sel = nodes$type == "T"
    )
  
  return(graph)
}

#' Measure diversity in a population
#'
#' @param population List of individuals
#' @return List of diversity metrics
#' @export
measure_diversity <- function(population) {
  # Number of unique phenotypes
  phenotypes <- sapply(population, function(ind) ind$phenotype)
  unique_phenotypes <- length(unique(phenotypes))
  
  # Mean genotype values and standard deviation for each non-terminal
  genotype_stats <- list()
  all_nts <- unique(unlist(lapply(population, function(ind) names(ind$genotype))))
  
  for (nt in all_nts) {
    # Extract values for this non-terminal
    values <- lapply(population, function(ind) {
      if (nt %in% names(ind$genotype)) ind$genotype[[nt]] else NULL
    })
    values <- unlist(values)
    
    if (length(values) > 0) {
      genotype_stats[[nt]] <- list(
        mean = mean(values),
        sd = sd(values),
        min = min(values),
        max = max(values)
      )
    }
  }
  
  # Phenotype length diversity
  lengths <- nchar(phenotypes)
  length_diversity <- list(
    mean = mean(lengths),
    sd = sd(lengths),
    min = min(lengths),
    max = max(lengths)
  )
  
  # Fitness diversity
  fitness <- sapply(population, function(ind) ind$fitness)
  fitness_diversity <- list(
    mean = mean(fitness),
    sd = sd(fitness),
    min = min(fitness),
    max = max(fitness)
  )
  
  return(list(
    population_size = length(population),
    unique_phenotypes = unique_phenotypes,
    phenotype_length = length_diversity,
    fitness = fitness_diversity,
    genotype_stats = genotype_stats
  ))
}

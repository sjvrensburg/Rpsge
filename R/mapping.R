#' Generate a random individual for PSGE
#'
#' @param grammar Grammar structure from read_pcfg_grammar()
#' @param max_depth Maximum depth for the generated tree
#' @return Genotype representation as a list of real values, one list per non-terminal
#' @export
generate_random_individual <- function(grammar, max_depth = 17) {
  # Initialize genotype with empty list for each non-terminal
  genotype <- list()
  for (nt in grammar$non_terminals) {
    genotype[[nt]] <- numeric(0)
  }

  # Start building the individual from the start symbol
  result <- recursive_individual_creation(
    grammar = grammar,
    genotype = genotype,
    symbol = grammar$start_symbol,
    current_depth = 0,
    max_depth = max_depth
  )

  return(result$genotype)
}

#' Safely select a rule based on codon value and rule probabilities
#'
#' @param rules List of rules with probabilities
#' @param codon Probabilistic value for selection (0-1)
#' @return Selected rule index
#' @keywords internal
safe_rule_selection <- function(rules, codon) {
  # Ensure codon is valid
  if (is.null(codon) || is.na(codon) || !is.numeric(codon)) {
    codon <- runif(1)
  }

  # Bound codon to [0,1]
  codon <- max(0, min(1, codon))

  # Get probabilities
  probs <- sapply(rules, function(r) r$prob)

  # Ensure probabilities sum to 1
  if (sum(probs) == 0) {
    # Equal probability if all are zero
    probs <- rep(1/length(probs), length(probs))
  } else {
    probs <- probs / sum(probs)
  }

  # Select rule
  cum_prob <- 0
  selected_index <- length(rules)  # Default to last rule

  for (i in seq_along(probs)) {
    cum_prob <- cum_prob + probs[i]
    if (codon <= cum_prob) {
      selected_index <- i
      break
    }
  }

  return(selected_index)
}

#' Recursively create a random individual
#'
#' @param grammar Grammar structure
#' @param genotype Current genotype being built
#' @param symbol Current symbol to expand
#' @param current_depth Current tree depth
#' @param max_depth Maximum allowed tree depth
#' @return Updated genotype and reached depth
#' @keywords internal
recursive_individual_creation <- function(
    grammar, genotype, symbol, current_depth, max_depth) {
  # Guard against exceeding max depth
  if (current_depth > max_depth) {
    return(list(genotype = genotype, depth = current_depth))
  }

  # If we've reached maximum depth, use non-recursive rules
  if (current_depth >= max_depth) {
    # Get non-recursive rules
    grammar_with_flags <- identify_recursive_rules(grammar)
    non_recursive_rules <- get_non_recursive_rules(grammar_with_flags, symbol)

    if (length(non_recursive_rules) == 0) {
      # Fall back to any rule at max depth
      non_recursive_rules <- grammar$rules[[symbol]]
    }

    # Generate random value
    codon <- runif(1)
    genotype[[symbol]] <- c(genotype[[symbol]], codon)

    # Select rule using safe selection
    selected_index <- safe_rule_selection(non_recursive_rules, codon)

    # Map selected index from non-recursive rules to the original rule index
    original_indices <- which(sapply(grammar$rules[[symbol]], function(r) {
      identical(r, non_recursive_rules[[selected_index]])
    }))

    if (length(original_indices) == 0) {
      # If no match found (shouldn't happen but safety first)
      selected_rule <- grammar$rules[[symbol]][[1]]
    } else {
      selected_rule <- grammar$rules[[symbol]][[original_indices[1]]]
    }
  } else {
    # Normal case: select any rule based on probabilities
    codon <- runif(1)
    genotype[[symbol]] <- c(genotype[[symbol]], codon)

    # Select rule using safe selection
    selected_index <- safe_rule_selection(grammar$rules[[symbol]], codon)
    selected_rule <- grammar$rules[[symbol]][[selected_index]]
  }

  # Process selected rule's symbols for further expansion
  max_depth_reached <- current_depth
  for (i in seq_along(selected_rule$symbols)) {
    if (selected_rule$types[i] == "NT") {
      next_symbol <- selected_rule$symbols[i]
      # Only recurse if we haven't exceeded max depth
      if (current_depth + 1 <= max_depth) {
        result <- recursive_individual_creation(
          grammar, genotype, next_symbol, current_depth + 1, max_depth
        )
        genotype <- result$genotype
        max_depth_reached <- max(max_depth_reached, result$depth)
      }
    }
  }

  return(list(genotype = genotype, depth = max_depth_reached))
}

#' Map a genotype to a phenotype using PSGE principles
#'
#' @param grammar Grammar structure from read_pcfg_grammar()
#' @param genotype Genotype representation (list of real values per non-terminal)
#' @param max_depth Maximum tree depth allowed
#' @param seed Optional random seed for reproducible mapping (default: NULL)
#' @return List containing phenotype string and mapping information
#' @export
map_genotype <- function(grammar, genotype, max_depth = 17, seed = NULL) {
  # Set fixed seed for tests if none provided (this ensures consistency for same genotype)
  if (is.null(seed) && any(sapply(genotype, length) > 0)) {
    # Create deterministic seed from genotype content
    genotype_str <- paste(unlist(genotype), collapse = ";")
    seed <- sum(utf8ToInt(genotype_str))
  }

  # Save original random state
  if (exists(".Random.seed")) {
    old_seed <- .Random.seed
    restore_seed <- TRUE
  } else {
    restore_seed <- FALSE
  }

  # Set seed if provided or calculated
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create a complete copy of the genotype
  genotype_copy <- list()
  for (nt in grammar$non_terminals) {
    if (nt %in% names(genotype) && length(genotype[[nt]]) > 0) {
      genotype_copy[[nt]] <- genotype[[nt]]
    } else {
      # Initialize with empty vector
      genotype_copy[[nt]] <- numeric(0)
    }
  }

  # Initialize mapping positions for each non-terminal
  positions <- list()
  for (nt in grammar$non_terminals) {
    positions[[nt]] <- 1
  }

  # Start mapping from the start symbol
  result <- recursive_mapping(
    grammar = grammar,
    genotype = genotype_copy,
    positions = positions,
    symbol = grammar$start_symbol,
    current_depth = 0,
    max_depth = max_depth,
    use_position_seeds = TRUE  # Always use position-based seeds for consistency
  )

  # Restore original random state if needed
  if (restore_seed) {
    .Random.seed <- old_seed
  }

  # Clean up any placeholder text
  clean_output <- gsub("max_depth_exceeded|at_max_depth", "", result$output)

  return(list(
    phenotype = paste(clean_output, collapse = ""),
    depth = result$depth,
    positions = result$positions
  ))
}

#' Recursively map a genotype to a phenotype
#'
#' @param grammar Grammar structure
#' @param genotype Genotype representation
#' @param positions Current position in each non-terminal's list
#' @param symbol Current symbol to expand
#' @param current_depth Current tree depth
#' @param max_depth Maximum allowed tree depth
#' @param output Current output being built
#' @param use_position_seeds Whether to use position-based seeds for consistency
#' @return Updated output, depth reached, positions, and genotype
#' @keywords internal
recursive_mapping <- function(grammar, genotype, positions, symbol,
                              current_depth, max_depth, output = character(0),
                              use_position_seeds = FALSE,
                              debug_output = FALSE) {
  # Add debug output before calling C++
  if (debug_output) {
    cat("\nR wrapper - Initial state:\n")
    cat("Positions:\n")
    str(positions)
    cat("Genotype:\n")
    str(genotype)
  }

  # Call the C++ implementation
  result <- recursiveMappingCpp(grammar, genotype, positions, symbol,
                                current_depth, max_depth, output,
                                use_position_seeds, debug_output)

  # Return all updated values including genotype
  return(list(
    output = result$output,
    depth = result$depth,
    positions = result$positions,
    genotype = result$genotype
  ))
}

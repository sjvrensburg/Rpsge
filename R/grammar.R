# We want to express grammars as simple lists of lists.
# I.e., avoid the convoluted class-based implementation
# from previous versions.

#' Parse and validate a BNF grammar file for PSGE
#'
#' @param grammar_file Path to the BNF grammar file
#' @param initial_prob_method Method to initialize probabilities ("uniform" or "specified")
#' @param prob_map Optional list of probability maps for each non-terminal
#' @param validate Boolean indicating whether to validate the grammar with `validate_grammar`
#' @return A validated PCFG grammar structure
#' @export
read_pcfg_grammar <- function(grammar_file, initial_prob_method = "uniform", prob_map = NULL, validate = TRUE) {
  if (!file.exists(grammar_file)) {
    stop("Grammar file not found: ", grammar_file)
  }

  # Read file content
  lines <- readLines(grammar_file)

  # Initialize grammar structure
  grammar <- list(
    non_terminals = character(0),
    terminals = character(0),
    start_symbol = NULL,
    rules = list()
  )

  # Regular expressions for parsing
  nt_pattern <- "<[^>]+>"
  rule_sep <- "::="
  prod_sep <- "\\|"

  # Process each line
  for (line in lines) {
    line <- trimws(line)
    if (line == "" || startsWith(line, "#")) {
      next # Skip empty lines and comments
    }

    # Split into LHS and RHS
    parts <- strsplit(line, rule_sep, fixed = FALSE)[[1]]
    if (length(parts) != 2) {
      next # Skip invalid lines
    }

    # Extract non-terminal (left-hand side)
    lhs <- trimws(parts[1])
    if (!grepl(paste0("^", nt_pattern, "$"), lhs)) {
      warning("Invalid non-terminal format: ", lhs)
      next
    }

    # Clean non-terminal name (remove < >)
    nt_name <- gsub("[<>]", "", lhs)

    # Add to non-terminals list if new
    if (!(nt_name %in% grammar$non_terminals)) {
      grammar$non_terminals <- c(grammar$non_terminals, nt_name)
      grammar$rules[[nt_name]] <- list()
    }

    # Set start symbol if not set
    if (is.null(grammar$start_symbol)) {
      grammar$start_symbol <- nt_name
    }

    # Process right-hand side (production rules)
    rhs <- trimws(parts[2])
    productions <- strsplit(rhs, prod_sep, fixed = FALSE)[[1]]
    productions <- trimws(productions)

    # Process each production
    for (prod in productions) {
      # Parse tokens
      tokens <- extract_tokens(prod)

      # Add production to grammar
      prod_entry <- list(
        symbols = tokens$symbols,
        types = tokens$types
      )

      grammar$rules[[nt_name]] <- c(grammar$rules[[nt_name]], list(prod_entry))

      # Update terminals set
      for (i in seq_along(tokens$symbols)) {
        if (tokens$types[i] == "T" && !(tokens$symbols[i] %in% grammar$terminals)) {
          grammar$terminals <- c(grammar$terminals, tokens$symbols[i])
        }
      }
    }
  }

  # Initialize probabilities
  grammar <- initialize_probabilities(grammar, method = initial_prob_method, prob_map = prob_map)

  # Validate the grammar if requested
  if (validate) {
    grammar <- validate_grammar(grammar)  # Store the returned grammar
  }

  return(grammar)
}

#' Extract tokens (terminals and non-terminals) from a production string
#'
#' @param production String containing a production rule
#' @return List with symbols and their types
#' @keywords internal
extract_tokens <- function(production) {
  nt_pattern <- "<[^>]*>"

  # Special case for empty non-terminal
  if (production == "<>") {
    return(list(symbols = "", types = "NT"))
  }

  # Find all non-terminals in the production
  nt_matches <- gregexpr(nt_pattern, production)
  nt_positions <- attr(nt_matches[[1]], "match.length")
  nt_starts <- as.vector(nt_matches[[1]])

  if (nt_starts[1] == -1) {
    # No non-terminals, this is a terminal production
    return(list(
      symbols = trimws(production),
      types = "T"
    ))
  }

  # Process mixed terminals and non-terminals
  tokens <- list()
  symbols <- character(0)
  types <- character(0)

  # Current position in the string
  pos <- 1
  for (i in seq_along(nt_starts)) {
    if (nt_starts[i] > pos) {
      # There's a terminal before this non-terminal
      term <- trimws(substr(production, pos, nt_starts[i] - 1))
      if (term != "") {
        symbols <- c(symbols, term)
        types <- c(types, "T")
      }
    }

    # Extract the non-terminal
    nt <- substr(production, nt_starts[i], nt_starts[i] + nt_positions[i] - 1)
    nt_clean <- gsub("[<>]", "", nt)
    symbols <- c(symbols, nt_clean)
    types <- c(types, "NT")

    # Update position
    pos <- nt_starts[i] + nt_positions[i]
  }

  # Check if there's a terminal at the end
  if (pos <= nchar(production)) {
    term <- trimws(substr(production, pos, nchar(production)))
    if (term != "") {
      symbols <- c(symbols, term)
      types <- c(types, "T")
    }
  }

  return(list(symbols = symbols, types = types))
}

#' Initialize probabilities for grammar rules
#'
#' @param grammar Grammar structure to initialize
#' @param method Method to use ("uniform" or "specified")
#' @param prob_map Map of probabilities (required if method is "specified")
#' @return Grammar with initialized probabilities
#' @export
initialize_probabilities <- function(grammar, method = "uniform", prob_map = NULL) {
  if (method == "uniform") {
    # Assign uniform probabilities to each production rule
    for (nt in names(grammar$rules)) {
      n_rules <- length(grammar$rules[[nt]])
      for (i in seq_len(n_rules)) {
        grammar$rules[[nt]][[i]]$prob <- 1/n_rules
      }
    }
  } else if (method == "specified") {
    if (is.null(prob_map)) {
      stop("prob_map must be provided when method is 'specified'")
    }

    # Use specified probabilities
    for (nt in names(prob_map)) {
      if (!(nt %in% names(grammar$rules))) {
        warning("Non-terminal '", nt, "' in prob_map not found in grammar")
        next
      }

      probs <- prob_map[[nt]]
      if (length(probs) != length(grammar$rules[[nt]])) {
        warning("Number of probabilities for '", nt, "' doesn't match number of rules")
        next
      }

      # Assign specified probabilities
      for (i in seq_along(probs)) {
        grammar$rules[[nt]][[i]]$prob <- probs[i]
      }
    }

    # For any non-terminals not in prob_map, use uniform probabilities
    for (nt in names(grammar$rules)) {
      if (!(nt %in% names(prob_map))) {
        n_rules <- length(grammar$rules[[nt]])
        for (i in seq_len(n_rules)) {
          grammar$rules[[nt]][[i]]$prob <- 1/n_rules
        }
      }
    }
  } else {
    stop("Unknown probability initialization method: ", method)
  }

  return(grammar)
}

#' Validate a grammar structure
#'
#' @param grammar Grammar structure to validate
#' @param tolerance Tolerance for probability sum (default: 1e-10)
#' @return The validated grammar object if valid, error if invalid
#' @export
validate_grammar <- function(grammar, tolerance = 1e-10) {
  # Check required components in exact order to match tests
  if (is.null(grammar$non_terminals) || length(grammar$non_terminals) == 0) {
    stop("Grammar has no non-terminals")
  }

  if (is.null(grammar$rules) || length(grammar$rules) == 0) {
    stop("Grammar has no rules")
  }

  if (is.null(grammar$start_symbol)) {
    stop("Grammar has no start symbol")
  }

  if (!(grammar$start_symbol %in% grammar$non_terminals)) {
    stop("Start symbol '", grammar$start_symbol, "' is not a defined non-terminal")
  }

  # Check each non-terminal has rules
  for (nt in grammar$non_terminals) {
    if (!(nt %in% names(grammar$rules)) || length(grammar$rules[[nt]]) == 0) {
      stop("Non-terminal '", nt, "' has no production rules")
    }

    for (i in seq_along(grammar$rules[[nt]])) {
      rule <- grammar$rules[[nt]][[i]]

      # 1. Check for empty symbols
      if (is.null(rule$symbols) || length(rule$symbols) == 0) {
        stop("Rule ", i, " for '", nt, "' has no symbols")
      }

      # 2. Check for mismatched lengths
      if (is.null(rule$types) || length(rule$types) != length(rule$symbols)) {
        stop("Rule ", i, " for '", nt, "' has mismatched symbols and types")
      }

      # 3. Check for invalid symbol types
      if (!all(rule$types %in% c("T", "NT"))) {
        stop("Rule ", i, " for '", nt, "' has invalid symbol types")
      }

      # 4. Check for invalid probability - only flag negative probabilities
      if (is.null(rule$prob) || rule$prob < 0) {
        stop("Rule ", i, " for '", nt, "' has invalid probability: ", rule$prob)
      }

      # 5. Verify all non-terminals in the rule exist in the grammar
      for (j in seq_along(rule$symbols)) {
        if (rule$types[j] == "NT" && !(rule$symbols[j] %in% grammar$non_terminals)) {
          stop("Rule ", i, " for '", nt, "' references undefined non-terminal: ", rule$symbols[j])
        }
      }
    }
  }

  # Now, check and normalize probability distributions
  for (nt in grammar$non_terminals) {
    # Validate probability distribution
    probs <- sapply(grammar$rules[[nt]], function(rule) rule$prob)
    prob_sum <- sum(probs)

    # Handle normalization threshold with very specific logic
    if (abs(prob_sum - 1) > tolerance) {
      # Instead of using <= comparison which can fail at exact boundaries,
      # use an open interval (< 0.1) plus a separate check for exact equality
      if (abs(prob_sum - 1) < 0.1 || isTRUE(all.equal(abs(prob_sum - 1), 0.1, tolerance=1e-14))) {
        warning("Probabilities for '", nt, "' sum to ", prob_sum, ". Normalizing.")

        # Normalize the probabilities - update the grammar object
        for (i in seq_along(grammar$rules[[nt]])) {
          grammar$rules[[nt]][[i]]$prob <- grammar$rules[[nt]][[i]]$prob / prob_sum
        }
      } else {
        stop("Probabilities for '", nt, "' sum to ", prob_sum, ", not 1")
      }
    }
  }

  # Check for unreachable non-terminals
  reachable <- check_reachability(grammar)
  unreachable <- setdiff(grammar$non_terminals, reachable)
  if (length(unreachable) > 0) {
    warning("Grammar has unreachable non-terminals: ", paste(unreachable, collapse=", "))
  }

  return(grammar)
}

#' Check which non-terminals are reachable from the start symbol
#'
#' @param grammar Grammar structure to check
#' @return Character vector of reachable non-terminals
#' @keywords internal
check_reachability <- function(grammar) {
  reachable <- grammar$start_symbol
  visited <- character(0)

  while (length(reachable) > 0) {
    current <- reachable[1]
    reachable <- reachable[-1]
    visited <- c(visited, current)

    # Find all rules for the current non-terminal
    if (current %in% names(grammar$rules)) {
      for (rule in grammar$rules[[current]]) {
        # Find all non-terminals in this rule
        for (i in seq_along(rule$symbols)) {
          if (rule$types[i] == "NT") {
            nt <- rule$symbols[i]
            if (!(nt %in% visited) && !(nt %in% reachable)) {
              reachable <- c(reachable, nt)
            }
          }
        }
      }
    }
  }

  return(visited)
}

#' Print a grammar in BNF format with probabilities
#'
#' @param grammar Grammar structure to print
#' @param include_probs Whether to include probabilities in output
#' @return Invisibly returns the grammar
#' @export
print_grammar <- function(grammar, include_probs = TRUE) {
  cat("Grammar with", length(grammar$non_terminals), "non-terminals and",
      length(grammar$terminals), "terminals\n")
  cat("Start symbol:", grammar$start_symbol, "\n\n")

  for (nt in grammar$non_terminals) {
    cat("<", nt, "> ::= ", sep="")
    rules <- grammar$rules[[nt]]

    for (i in seq_along(rules)) {
      if (i > 1) cat(" | ")

      # Print symbols
      for (j in seq_along(rules[[i]]$symbols)) {
        if (rules[[i]]$types[j] == "NT") {
          cat("<", rules[[i]]$symbols[j], ">", sep="")
        } else {
          cat(rules[[i]]$symbols[j])
        }
        if (j < length(rules[[i]]$symbols)) cat(" ")
      }

      # Print probability if requested
      if (include_probs) {
        cat(" [", formatC(rules[[i]]$prob, digits=4, format="f"), "]", sep="")
      }
    }
    cat("\n")
  }

  invisible(grammar)
}

#' Get rules for a specific non-terminal
#'
#' @param grammar Grammar structure
#' @param non_terminal Name of the non-terminal
#' @return List of rules for the non-terminal
#' @export
get_rules <- function(grammar, non_terminal) {
  if (!(non_terminal %in% names(grammar$rules))) {
    stop("Non-terminal '", non_terminal, "' not found in grammar")
  }
  return(grammar$rules[[non_terminal]])
}

#' Set probabilities for rules of a non-terminal
#'
#' @param grammar Grammar structure
#' @param non_terminal Name of the non-terminal
#' @param probabilities Vector of probabilities
#' @return Updated grammar
#' @export
set_probabilities <- function(grammar, non_terminal, probabilities) {
  if (!(non_terminal %in% names(grammar$rules))) {
    stop("Non-terminal '", non_terminal, "' not found in grammar")
  }

  rules <- grammar$rules[[non_terminal]]
  if (length(probabilities) != length(rules)) {
    stop("Number of probabilities (", length(probabilities),
         ") doesn't match number of rules (", length(rules), ") for '", non_terminal, "'")
  }

  if (abs(sum(probabilities) - 1) > 1e-10) {
    stop("Probabilities must sum to 1, got ", sum(probabilities))
  }

  for (i in seq_along(rules)) {
    grammar$rules[[non_terminal]][[i]]$prob <- probabilities[i]
  }

  return(grammar)
}

#' Add a new rule to a non-terminal
#'
#' @param grammar Grammar structure
#' @param non_terminal Name of the non-terminal
#' @param rule_str String representation of the rule
#' @param prob Probability for the new rule (if NULL, probabilities will be redistributed)
#' @return Updated grammar
#' @export
add_rule <- function(grammar, non_terminal, rule_str, prob = NULL) {
  # Check if non-terminal exists
  if (!(non_terminal %in% grammar$non_terminals)) {
    # Add new non-terminal
    grammar$non_terminals <- c(grammar$non_terminals, non_terminal)
    grammar$rules[[non_terminal]] <- list()
  }

  # Parse the rule
  tokens <- extract_tokens(rule_str)
  new_rule <- list(
    symbols = tokens$symbols,
    types = tokens$types
  )

  # Calculate probabilities
  current_rules <- grammar$rules[[non_terminal]]
  n_rules <- length(current_rules)

  if (is.null(prob)) {
    # Distribute probabilities uniformly
    new_prob <- 1 / (n_rules + 1)
    for (i in seq_len(n_rules)) {
      grammar$rules[[non_terminal]][[i]]$prob <- new_prob
    }
    new_rule$prob <- new_prob
  } else {
    if (prob <= 0 || prob >= 1) {
      stop("Probability must be between 0 and 1")
    }

    # Redistribute remaining probability
    remaining <- 1 - prob
    if (n_rules > 0) {
      current_total <- sum(sapply(current_rules, function(r) r$prob))
      scale_factor <- remaining / current_total
      for (i in seq_len(n_rules)) {
        grammar$rules[[non_terminal]][[i]]$prob <- grammar$rules[[non_terminal]][[i]]$prob * scale_factor
      }
    }
    new_rule$prob <- prob
  }

  # Add the new rule
  grammar$rules[[non_terminal]] <- c(grammar$rules[[non_terminal]], list(new_rule))

  # Update terminals if needed
  for (i in seq_along(tokens$symbols)) {
    if (tokens$types[i] == "T" && !(tokens$symbols[i] %in% grammar$terminals)) {
      grammar$terminals <- c(grammar$terminals, tokens$symbols[i])
    }
  }

  # Validate the modified grammar
  grammar <- validate_grammar(grammar)  # Store the returned grammar

  return(grammar)
}

#' Create an empty grammar structure
#'
#' @param start Name of the start symbol
#' @return Empty grammar structure
#' @export
new_grammar <- function(start) {
  grammar <- list(
    non_terminals = start,  # Add start symbol to non-terminals
    terminals = character(0),
    start_symbol = start,
    rules = list()
  )
  # Initialize empty rules list for the start symbol
  grammar$rules[[start]] <- list()
  return(grammar)
}

#' Identify recursive and non-recursive rules
#'
#' @param grammar Grammar structure
#' @return Updated grammar with recursive_rules flags
#' @export
identify_recursive_rules <- function(grammar) {
  for (nt in names(grammar$rules)) {
    for (i in seq_along(grammar$rules[[nt]])) {
      rule <- grammar$rules[[nt]][[i]]
      is_recursive <- FALSE

      for (j in seq_along(rule$symbols)) {
        if (rule$types[j] == "NT" && rule$symbols[j] == nt) {
          is_recursive <- TRUE
          break
        }
      }

      grammar$rules[[nt]][[i]]$is_recursive <- is_recursive
    }
  }

  return(grammar)
}

#' Get non-recursive rules for a non-terminal
#'
#' @param grammar Grammar structure
#' @param non_terminal Name of the non-terminal
#' @return List of non-recursive rules
#' @export
get_non_recursive_rules <- function(grammar, non_terminal) {
  if (!(non_terminal %in% names(grammar$rules))) {
    stop("Non-terminal '", non_terminal, "' not found in grammar")
  }

  # Ensure recursive flags are set
  if (!all(sapply(grammar$rules[[non_terminal]], function(r) !is.null(r$is_recursive)))) {
    grammar <- identify_recursive_rules(grammar)
  }

  # Filter non-recursive rules
  non_recursive <- list()
  for (rule in grammar$rules[[non_terminal]]) {
    if (!rule$is_recursive) {
      non_recursive <- c(non_recursive, list(rule))
    }
  }

  if (length(non_recursive) == 0) {
    warning("No non-recursive rules found for '", non_terminal, "'")
  }

  return(non_recursive)
}

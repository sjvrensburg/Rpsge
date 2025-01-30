#' Grammar Class for PSGE
#'
#' @description
#' An R6 class implementing a probabilistic context-free grammar (PCFG) for use
#' with Probabilistic Structured Grammatical Evolution.
#'
#' @details
#' This class handles reading and parsing BNF grammars, managing production rules,
#' and updating rule probabilities during evolution. It supports both standard
#' context-free grammars and probabilistic context-free grammars.
#'
#' @examples
#' \dontrun{
#' # Create grammar from file
#' grammar <- Grammar$new("path/to/grammar.bnf")
#'
#' # Generate uniform probabilities
#' grammar$generate_uniform_pcfg()
#'
#' # Update probabilities based on an individual
#' grammar$update_pcfg_probabilities(individual, learning_factor = 0.01)
#' }
#'
#' @export
Grammar <- R6::R6Class(
  "Grammar",

  private = list(
    # Constants as private fields
    .NT = "NT",
    .T = "T",
    .NT_PATTERN = "<.+?>",
    .RULE_SEPARATOR = "::=",
    .PRODUCTION_SEPARATOR = "\\|",

    # Private fields
    .grammar_file = NULL,
    .grammar = NULL,
    .non_terminals = NULL,
    .terminals = NULL,
    .ordered_non_terminals = NULL,
    .start_rule = NULL,
    .max_depth = NULL,
    .max_init_depth = NULL,
    .pcfg = NULL,
    .pcfg_mask = NULL,
    .pcfg_path = NULL,
    .index_of_non_terminal = NULL,
    .shortest_path = NULL,

    # Cleanup method
    .finalize = function() {
      private$.grammar <- NULL
      private$.pcfg <- NULL
      private$.shortest_path <- NULL
      gc() # Suggest garbage collection
    },

    # Private helper functions
    .validate_grammar = function() {
      # Validate grammar structure
      if (length(private$.non_terminals) == 0) {
        stop("Grammar must have at least one non-terminal")
      }
      if (is.null(private$.start_rule)) {
        stop("Grammar must have a start rule")
      }
      # Check for unreachable rules
      reachable <- self$find_reachable_rules()
      unreachable <- setdiff(names(private$.grammar), reachable)
      if (length(unreachable) > 0) {
        warning("Grammar contains unreachable rules: ",
                paste(unreachable, collapse = ", "))
      }
    },

    #' @description Recursive mapping helper
    #' @param mapping_rules List of production rules
    #' @param positions_to_map Vector tracking positions
    #' @param current_sym Current symbol being processed
    #' @param current_depth Current tree depth
    #' @param output Character vector collecting output
    #' @return Maximum depth reached
    .recursive_mapping = function(mapping_rules, positions_to_map,
                                  current_sym, current_depth, output) {
      if (is.null(private$.max_depth)) {
        private$.max_depth <- 17  # Default value
        warning("max_depth was not initialized, using default value of 17")
      }

      depths <- current_depth

      if (current_sym[[2]] == self$T) {
        # Terminal symbol - add to output
        output$symbols <- c(output$symbols, current_sym[[1]])
      } else {
        # Non-terminal processing
        current_sym_pos <- which(private$.non_terminals == current_sym[[1]])
        choices <- private$.grammar[[current_sym[[1]]]]
        nt_index <- private$.index_of_non_terminal[[current_sym[[1]]]]

        # Generate or use existing codon
        if (positions_to_map[current_sym_pos] >=
            length(mapping_rules[[current_sym_pos]])) {

          # Generate new codon
          codon <- stats::runif(1)

          # Handle depth limit case
          if (current_depth > private$.max_depth) {
            shortest_rules <- private$.shortest_path[[current_sym[[1]]]]$rules
            expansion_possibility <- private$.select_non_recursive_rule(
              shortest_rules, codon, nt_index
            )
          } else {
            # Normal expansion
            expansion_possibility <- private$.select_rule(
              codon, nt_index, length(choices)
            )
          }

          # Add new mapping rule
          mapping_rules[[current_sym_pos]] <- c(
            mapping_rules[[current_sym_pos]],
            list(c(expansion_possibility - 1, codon, current_depth))
          )
        } else {
          # Use existing codon but recalculate rule based on current probabilities
          codon <- mapping_rules[[current_sym_pos]][[
            positions_to_map[current_sym_pos] + 1
          ]][2]

          if (current_depth > private$.max_depth) {
            shortest_rules <- private$.shortest_path[[current_sym[[1]]]]$rules
            expansion_possibility <- private$.select_non_recursive_rule(
              shortest_rules, codon, nt_index
            )
          } else {
            expansion_possibility <- private$.select_rule(
              codon, nt_index, length(choices)
            )
          }

          # Update existing rule
          mapping_rules[[current_sym_pos]][[
            positions_to_map[current_sym_pos] + 1
          ]] <- c(expansion_possibility - 1, codon, current_depth)
        }

        # Update position counter
        positions_to_map[current_sym_pos] <- positions_to_map[current_sym_pos] + 1

        # Process chosen production
        next_to_expand <- choices[[expansion_possibility]]
        for (next_sym in next_to_expand) {
          sub_depth <- private$.recursive_mapping(
            mapping_rules = mapping_rules,
            positions_to_map = positions_to_map,
            current_sym = next_sym,
            current_depth = current_depth + 1,
            output = output
          )
          depths <- c(depths, sub_depth)
        }
      }

      return(max(depths))
    },

    #' @description Select production rule based on probabilities
    #' @param codon Random value to use for selection
    #' @param nt_index Non-terminal index
    #' @param n_choices Number of choices available
    #' @return Selected rule index
    .select_rule = function(codon, nt_index, n_choices) {
      prob_sum <- 0
      for (i in seq_len(n_choices)) {
        prob_sum <- prob_sum + private$.pcfg[nt_index, i]
        if (codon <= prob_sum) {
          return(i)
        }
      }
      return(n_choices)  # Default to last rule if no match
    },

    #' @description Select non-recursive rule based on probabilities
    #' @param shortest_rules List of valid non-recursive rules
    #' @param codon Random value to use for selection
    #' @param nt_index Non-terminal index
    #' @return Selected rule index
    .select_non_recursive_rule = function(shortest_rules, codon, nt_index) {
      # Validate input
      if (length(shortest_rules) == 0) {
        stop("No non-recursive rules provided")
      }

      # Extract rules in correct format
      rules <- shortest_rules[[1]]  # Get first (and should be only) rule set
      if (!is.list(rules) || !("index" %in% names(rules))) {
        stop("Invalid rule format in shortest_rules")
      }

      # Get rule index
      return(rules$index)
    },

    #' @description Apply Python-specific filtering
    #' @param text Text to filter
    #' @return Filtered text
    .python_filter = function(text) {
      # Replace operators
      text <- gsub("\\\\le", "<=", text)
      text <- gsub("\\\\ge", ">=", text)
      text <- gsub("\\\\l", "<", text)
      text <- gsub("\\\\g", ">", text)
      text <- gsub("\\\\eb", "|", text)

      return(text)
    }
  ),

  public = list(
    #' @description
    #' Create a new Grammar instance
    #'
    #' @param grammar_path Character path to grammar file (optional)
    #' @return A new `Grammar` object
    initialize = function(grammar_path = NULL) {
      private$.max_init_depth <- 6
      private$.max_depth <- 17      # Default value matching Python implementation
      private$.grammar_file <- grammar_path
      private$.grammar <- list()
      private$.non_terminals <- character()
      private$.terminals <- character()
      private$.start_rule <- NULL
      private$.pcfg <- NULL
      private$.pcfg_mask <- NULL

      if (!is.null(grammar_path)) {
        self$read_grammar()
      }
    },

    #' @description
    #' Read grammar from file
    #'
    #' @return Invisible self (for method chaining)
    read_grammar = function() {
      if (is.null(private$.grammar_file)) {
        stop("Grammar file path not specified")
      }

      # Read file safely
      lines <- tryCatch({
        readLines(private$.grammar_file)
      }, error = function(e) {
        stop("Failed to read grammar file: ", e$message)
      })

      # Filter comments and empty lines
      lines <- lines[!grepl("^#", lines) &
                       nzchar(trimws(lines))]

      # Parse rules
      for (line in lines) {
        if (grepl("::=", line)) {
          parts <- strsplit(line, "::=")[[1]]
          left_side <- trimws(parts[1])

          # Validate left side
          if (!grepl("<.+?>", left_side)) {
            stop("Invalid non-terminal format: ", left_side)
          }

          private$.non_terminals <- c(private$.non_terminals, left_side)

          # Set start rule if not set
          if (is.null(private$.start_rule)) {
            private$.start_rule <- c(left_side, self$NT)
          }

          # Parse productions more efficiently
          productions <- strsplit(parts[2], private$.PRODUCTION_SEPARATOR)[[1]]
          temp_productions <- lapply(trimws(productions), function(prod) {
            if (!grepl(private$.NT_PATTERN, prod)) {
              # Simple terminal production
              if (prod == "None") prod <- ""
              private$.terminals <- c(private$.terminals, prod)
              list(list(prod, private$.T))  # Single symbol production
            } else {
              # Complex production with multiple symbols
              symbols <- strsplit(prod, "\\s+")[[1]]
              symbols <- symbols[nzchar(trimws(symbols))]  # Remove empty strings

              # Process each symbol
              lapply(symbols, function(val) {
                val <- trimws(val)
                if (!grepl(private$.NT_PATTERN, val)) {
                  private$.terminals <- c(private$.terminals, val)
                  list(val, private$.T)
                } else {
                  list(val, private$.NT)
                }
              })
            }
          })

          private$.grammar[[left_side]] <- temp_productions
        }
      }

      # Remove duplicates
      private$.terminals <- unique(private$.terminals)
      private$.non_terminals <- unique(private$.non_terminals)
      private$.ordered_non_terminals <- private$.non_terminals  # Maintain order

      # Update index mapping
      private$.index_of_non_terminal <- stats::setNames(
        seq_along(private$.non_terminals),
        private$.non_terminals
      )

      # Validate and initialize
      private$.validate_grammar()
      self$generate_uniform_pcfg()
      self$find_shortest_path()
    },

    #' @description
    #' Save grammar state to file
    #'
    #' @param file_path Character path where to save grammar
    #'
    #' @return Invisible self (for method chaining)
    save_grammar = function(file_path) {
      grammar_data <- list(
        grammar = private$.grammar,
        non_terminals = private$.non_terminals,
        terminals = private$.terminals,
        pcfg = private$.pcfg,
        pcfg_mask = private$.pcfg_mask,
        shortest_path = private$.shortest_path
      )
      saveRDS(grammar_data, file = file_path)
      invisible(self)
    },

    #' @description
    #' Load previously saved grammar state from file
    #'
    #' @param file_path Character path to the saved grammar file (.rds)
    #' @return Invisible self (for method chaining)
    #' @examples
    #' \dontrun{
    #' grammar <- Grammar$new()
    #' grammar$load_grammar("saved_grammar.rds")
    #' }
    load_grammar = function(file_path) {
      grammar_data <- readRDS(file = file_path)
      private$.grammar <- grammar_data$grammar
      private$.non_terminals <- grammar_data$non_terminals
      private$.terminals <- grammar_data$terminals
      private$.pcfg <- grammar_data$pcfg
      private$.pcfg_mask <- grammar_data$pcfg_mask
      private$.shortest_path <- grammar_data$shortest_path
      private$.index_of_non_terminal <- stats::setNames(
        seq_along(private$.non_terminals),
        private$.non_terminals
      )
      invisible(self)
    },

    #' @description Find reachable rules from start symbol
    find_reachable_rules = function() {
      reachable <- character()
      to_process <- private$.start_rule[1]

      while (length(to_process) > 0) {
        current <- to_process[1]
        to_process <- to_process[-1]

        if (!(current %in% reachable)) {
          reachable <- c(reachable, current)

          # Add all non-terminals from productions
          if (current %in% names(private$.grammar)) {
            for (prod in private$.grammar[[current]]) {
              for (sym in prod) {
                if (sym[[2]] == self$NT) {
                  to_process <- c(to_process, sym[[1]])
                }
              }
            }
          }
        }
      }
      reachable
    },

    #' @description Generate uniform probabilities for PCFG
    generate_uniform_pcfg = function() {
      # Get maximum number of production rules for any non-terminal
      max_rules <- max(vapply(private$.grammar, length, numeric(1)))
      n_non_terminals <- length(private$.non_terminals)

      # Create probability matrix
      private$.pcfg <- matrix(0, nrow = n_non_terminals, ncol = max_rules)
      private$.pcfg_mask <- matrix(FALSE, nrow = n_non_terminals, ncol = max_rules)

      # Generate index mapping for non-terminals if not exists
      if (length(private$.index_of_non_terminal) == 0) {
        private$.index_of_non_terminal <- stats::setNames(seq_along(private$.non_terminals),
                                                          private$.non_terminals)
      }

      # Fill probability matrix with uniform distributions
      for (i in seq_len(n_non_terminals)) {
        nt <- private$.non_terminals[i]
        n_prods <- length(private$.grammar[[nt]])
        prob <- 1 / n_prods

        # Set probabilities and mask
        private$.pcfg[i, seq_len(n_prods)] <- prob
        private$.pcfg_mask[i, seq_len(n_prods)] <- TRUE
      }

      invisible(self)
    },

    #' @description Update PCFG probabilities based on an individual
    #' @param individual List containing genotype and mapping information
    #' @param learning_factor Numeric learning rate for probability updates
    update_pcfg_probabilities = function(individual, learning_factor) {
      # Get expansion counts from individual
      grammar_counter <- self$get_grammar_counter(individual)

      # Update probabilities using C++ function
      private$.pcfg <- update_pcfg_probabilities(
        pcfg = private$.pcfg,
        mask = private$.pcfg_mask,
        grammar_counter = grammar_counter,
        learning_factor = learning_factor
      )

      invisible(self)
    },

    #' @description Count grammar rule usage in an individual
    #' @param individual List containing genotype and mapping information
    #' @return List of numeric vectors containing rule counts
    get_grammar_counter = function(individual) {
      # Initialize counters for each non-terminal
      n_nts <- length(private$.non_terminals)
      counter <- vector("list", n_nts)
      names(counter) <- private$.non_terminals

      # Set counter lengths based on number of rules
      for (i in seq_len(n_nts)) {
        nt <- private$.non_terminals[i]
        counter[[i]] <- numeric(length(private$.grammar[[nt]]))
      }

      # Count expansions from genotype up to mapping values
      for (i in seq_len(n_nts)) {
        nt <- private$.non_terminals[i]
        nt_idx <- private$.index_of_non_terminal[[nt]]
        expansions <- individual$genotype[[nt_idx]]

        # Only count up to used mapping values
        if (length(expansions) > 0) {
          mapped_count <- individual$mapping_values[nt_idx]
          if (mapped_count > 0) {
            # Extract rule indices and count them
            rule_indices <- vapply(expansions[seq_len(mapped_count)], function(x)
              x[1] + 1, numeric(1))
            tab <- table(factor(rule_indices, levels = seq_along(counter[[i]])))
            counter[[i]] <- as.numeric(tab)
          }
        }
      }

      return(counter)
    },

    #' @description Get probability for a specific rule
    #' @param nt_index Index of non-terminal
    #' @param rule_index Index of production rule
    #' @return Numeric probability value
    get_probability = function(nt_index, rule_index) {
      if (private$.pcfg_mask[nt_index, rule_index]) {
        return(private$.pcfg[nt_index, rule_index])
      }
      return(0)
    },

    #' @description Get all probabilities for a non-terminal
    #' @param nt_index Index of non-terminal
    #' @return Numeric vector of probabilities
    get_probabilities_non_terminal = function(nt_index) {
      return(private$.pcfg[nt_index, private$.pcfg_mask[nt_index, ]])
    },

    #' @description Find shortest paths for all non-terminals
    #' @return Invisible self
    find_shortest_path = function() {
      # Initialize shortest paths storage
      private$.shortest_path <- list()

      # Process each non-terminal
      for (nt in private$.non_terminals) {
        # Find shortest path for this non-terminal
        depth <- self$find_minimum_path(current_symbol = list(nt, self$NT),
                                        open_symbols = character())

        # Store result if valid path found
        if (!is.null(depth)) {
          private$.shortest_path[[nt]] <- depth
        }
      }

      invisible(self)
    },

    #' @description Calculate minimum path for a symbol
    #' @param current_symbol List containing symbol and type
    #' @param open_symbols Character vector of symbols being processed
    #' @return List containing path depth and rules, or NULL if no path found
    find_minimum_path = function(current_symbol, open_symbols) {
      # Terminal symbols have depth 0
      if (current_symbol[[2]] == self$T) {
        return(list(depth = 0, rules = list()))
      }

      # Track symbol to detect cycles
      symbol_name <- current_symbol[[1]]
      if (symbol_name %in% open_symbols) {
        return(NULL)  # Cycle detected
      }
      open_symbols <- c(open_symbols, symbol_name)

      best_depth <- Inf
      best_rules <- list()

      # Check each production rule
      for (rule_idx in seq_along(private$.grammar[[symbol_name]])) {
        rule <- private$.grammar[[symbol_name]][[rule_idx]]

        # Skip recursive rules (containing current non-terminal)
        is_recursive <- any(vapply(rule, function(sym) {
          sym[[1]] == symbol_name && sym[[2]] == self$NT
        }, logical(1)))

        if (!is_recursive) {
          # Calculate max depth of this rule
          max_rule_depth <- 0
          rule_valid <- TRUE

          # Process each symbol in rule
          for (sym in rule) {
            # Recursively find path for this symbol
            sub_path <- self$find_minimum_path(sym, open_symbols)

            if (is.null(sub_path)) {
              rule_valid <- FALSE
              break
            }

            max_rule_depth <- max(max_rule_depth, sub_path$depth)
          }

          # Update best path if this rule is better
          if (rule_valid && max_rule_depth < best_depth) {
            best_depth <- max_rule_depth
            best_rules <- list(
              rule = rule,
              index = rule_idx
            )
          }
        }
      }

      # Return NULL if no valid path found
      if (best_depth == Inf) {
        return(NULL)
      }

      # Return path info
      return(list(
        depth = best_depth + 1,
        rules = list(best_rules)  # Wrap in list to match expected structure
      ))
    },

    #' @description Get non-recursive production rules for a symbol
    #' @param symbol Character non-terminal symbol
    #' @return List containing rules and their total probability
    get_non_recursive_rules = function(symbol) {
      rules <- private$.grammar[[symbol]]
      nt_idx <- private$.index_of_non_terminal[[symbol]]

      # Find non-recursive rules
      non_recursive <- lapply(seq_along(rules), function(i) {
        rule <- rules[[i]]
        is_recursive <- any(vapply(rule, function(sym) {
          sym[[1]] == symbol && sym[[2]] == self$NT
        }, logical(1)))

        if (!is_recursive) {
          return(list(
            rule = rule,
            index = i,
            probability = private$.pcfg[nt_idx, i]
          ))
        }
        return(NULL)
      })

      # Filter out NULL entries and calculate total probability
      non_recursive <- Filter(Negate(is.null), non_recursive)
      total_prob <- sum(vapply(non_recursive, function(x) x$probability, numeric(1)))

      return(list(
        rules = non_recursive,
        total_probability = total_prob
      ))
    },

    #' @description Validate maximum depths for all paths
    #' @param max_depth Maximum allowed depth
    #' @return Logical indicating if all paths are valid
    validate_depths = function(max_depth) {
      all_valid <- TRUE

      for (nt in private$.non_terminals) {
        path <- private$.shortest_path[[nt]]
        if (!is.null(path) && path$depth > max_depth) {
          warning("Path for ", nt, " exceeds max depth: ", path$depth)
          all_valid <- FALSE
        }
      }

      return(all_valid)
    },

    #' @description Map genotype to phenotype
    #' @param mapping_rules List of production rules to apply
    #' @param positions_to_map Vector tracking mapping positions
    #' @param needs_python_filter Logical indicating if Python filtering needed
    #' @return List containing:
    #'   - Character string of mapped expression
    #'   - Numeric maximum depth reached
    mapping = function(mapping_rules, positions_to_map = NULL, needs_python_filter = FALSE) {
      # Initialize position tracking if not provided
      if (is.null(positions_to_map)) {
        positions_to_map <- numeric(length(private$.non_terminals))
      }

      # Initialize output collection as a reference object (environment)
      output <- new.env()
      output$symbols <- character()

      # Perform recursive mapping
      max_depth <- private$.recursive_mapping(
        mapping_rules = mapping_rules,
        positions_to_map = positions_to_map,
        current_sym = private$.start_rule,
        current_depth = 0,
        output = output
      )

      # Combine output into single string
      phenotype <- paste(output$symbols, collapse = "")

      # Apply Python filtering if needed
      if (needs_python_filter) {
        phenotype <- private$.python_filter(phenotype)
      }

      return(list(
        phenotype = phenotype,
        max_depth = max_depth
      ))
    },

    #' @description Create random valid individual
    #' @param genome List to store genotype
    #' @param symbol Starting symbol
    #' @param current_depth Current tree depth
    #' @return Maximum depth reached
    recursive_individual_creation = function(genome, symbol, current_depth) {
      # Add validation at the start
      if (is.null(private$.max_init_depth)) {
        private$.max_init_depth <- 6  # Default value
        warning("max_init_depth was not initialized, using default value of 6")
      }

      if (is.null(current_depth)) {
        current_depth <- 0
        warning("current_depth was NULL, defaulting to 0")
      }

      # Generate random codon
      codon <- stats::runif(1)
      nt_index <- private$.index_of_non_terminal[[symbol]]

      # Handle depth limit case
      if (current_depth > private$.max_init_depth) {
        # Get non-recursive options
        shortest_path <- private$.shortest_path[[symbol]]$rules
        if (length(shortest_path) == 0) {
          stop("No valid non-recursive path found for ", symbol)
        }

        # Calculate probabilities for non-recursive rules
        prob_non_recursive <- 0
        non_recursive_indices <- integer()

        for (rule in shortest_path) {
          rule_idx <- rule$index
          prob_non_recursive <- prob_non_recursive +
            private$.pcfg[nt_index, rule_idx]
          non_recursive_indices <- c(non_recursive_indices, rule_idx)
        }

        # Select rule based on normalized probabilities
        prob_sum <- 0
        expansion_possibility <- NULL

        for (rule_idx in non_recursive_indices) {
          new_prob <- private$.pcfg[nt_index, rule_idx] / prob_non_recursive
          prob_sum <- prob_sum + new_prob
          if (codon <= prob_sum) {
            expansion_possibility <- rule_idx
            break
          }
        }
      } else {
        # Normal expansion case - use full PCFG
        prob_sum <- 0
        for (rule_idx in seq_along(private$.grammar[[symbol]])) {
          prob_sum <- prob_sum + private$.pcfg[nt_index, rule_idx]
          if (codon <= prob_sum) {
            expansion_possibility <- rule_idx
            break
          }
        }
      }

      # Add expansion to genome
      nt_position <- which(private$.non_terminals == symbol)
      genome[[nt_position]] <- c(
        genome[[nt_position]],
        list(c(expansion_possibility - 1, codon, current_depth))
      )

      # Get expansion symbols
      expansion_symbols <- private$.grammar[[symbol]][[expansion_possibility]]

      # Track depths
      depths <- current_depth

      # Recurse for non-terminals
      for (sym in expansion_symbols) {
        if (sym[[2]] == self$NT) {
          depths <- c(depths,
                      self$recursive_individual_creation(
                        genome = genome,
                        symbol = sym[[1]],
                        current_depth = current_depth + 1
                      )
          )
        }
      }

      return(max(depths))
    }
  ),

  active = list(

    #' @field NT Get NT constant
    NT = function() private$.NT,

    #' @field T Get T constant
    T = function() private$.T,

    #' @field NT_PATTERN Get NT_PATTERN constant
    NT_PATTERN = function() private$.NT_PATTERN,

    #' @field RULE_SEPARATOR Get RULE_SEPARATOR constant
    RULE_SEPARATOR = function() private$.RULE_SEPARATOR,

    #' @field PRODUCTION_SEPARATOR Get PRODUCTION_SEPARATOR constant
    PRODUCTION_SEPARATOR = function() private$.PRODUCTION_SEPARATOR,

    #' @field grammar Get the grammar rules
    grammar = function() private$.grammar,

    #' @field non_terminals Get non-terminal symbols
    non_terminals = function() private$.non_terminals,

    #' @field terminals Get terminal symbols
    terminals = function() private$.terminals,

    #' @field pcfg Get probabilistic context-free grammar matrix
    pcfg = function() private$.pcfg
  )
)

#' Create a new Grammar instance
#'
#' @param grammar_path Character path to grammar file
#' @return A new Grammar object
#' @export
create_grammar <- function(grammar_path = NULL) {
  Grammar$new(grammar_path)
}

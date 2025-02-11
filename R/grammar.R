round_prob <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0) return(x)
  if (any(is.na(x))) return(x)  # Preserve NA values
  return(round(x, digits))
}

# R equivalent for Python's static methods (defined outside the R6 class)
r_filter <- function(txt, needs_r_filter) {
  # String replacements - still relevant
  txt <- gsub(r"(\\le)", "<=", txt)
  txt <- gsub(r"(\\ge)", ">=", txt)
  txt <- gsub(r"(\\l)", "<", txt)
  txt <- gsub(r"(\\g)", ">", txt)
  txt <- gsub(r"(\\eb)", "\\|", txt)
  txt <- gsub(r"(\\n)", "\n", txt)
  txt <- gsub(r"(\\t)", "\t", txt)
  txt <- gsub(r"(\{:)", "{", txt)
  txt <- gsub(r"(:\})", "}", txt)

  if (needs_r_filter) { # Apply formatting only when needed
    # **Remove markers BEFORE styling**
    txt_no_markers <- gsub(r"(\{:)", "", txt)
    txt_no_markers <- gsub(r"(:\})", "", txt_no_markers)
    txt_no_markers <- gsub(r"(\\n)", "\n", txt_no_markers)

    tryCatch({
      # Attempt to style the text as R code using styler::style_text
      formatted_txt <- styler::style_text(txt_no_markers)
      # styler::style_text returns a list of formatted code lines.
      txt <- paste(formatted_txt, collapse = "\n") # Convert to single string
    }, error = function(e) {
      # If styler fails to parse or style (e.g., if 'txt' is not valid R-like code)
      warning("styler formatting failed: ", e$message, ". Returning unstyled text (markers removed).")
      # In case of error, return the original text without styling
      return(txt_no_markers)
    })
  }
  return(trimws(txt))
}

Grammar <- R6::R6Class("Grammar",
                   public = list(
                     NT = "NT", # Class attribute NT = "NT"
                     T = "T",   # Class attribute T = "T"
                     NT_PATTERN = "(<.+?>)", # Class attribute NT_PATTERN = "(<.+?>)"
                     RULE_SEPARATOR = "::=", # Class attribute RULE_SEPARATOR = "::="
                     PRODUCTION_SEPARATOR = "|", # Class attribute PRODUCTION_SEPARATOR = "|"

                     # Initialize all class attributes
                     initialize = function() {
                       # Basic attributes
                       self$grammar_file <- NULL  # Path to grammar file
                       self$grammar <- list()     # Main grammar structure (Python dict equivalent)
                       self$productions_labels <- list()  # Production labels dictionary

                       # Sets for terminals and non-terminals
                       self$non_terminals <- OrderedUniqueSet$new()
                       self$terminals <- OrderedUniqueSet$new()
                       self$ordered_non_terminals <- OrderedUniqueSet$new()  # Ordered set of non-terminals

                       # Options and rules
                       self$non_recursive_options <- list()  # Non-recursive production options
                       self$number_of_options_by_non_terminal <- NULL  # Count of options per non-terminal
                       self$start_rule <- NULL  # Starting rule for the grammar

                       # Depth controls
                       self$max_depth <- NULL  # Maximum tree depth
                       self$max_init_depth <- NULL  # Maximum initial tree depth
                       self$max_number_prod_rules <- 0  # Maximum number of production rules

                       # PCFG related attributes
                       self$pcfg <- NULL  # Probabilistic Context-Free Grammar matrix
                       self$pcfg_mask <- NULL  # Mask for PCFG
                       self$pcfg_path <- NULL  # Path to PCFG file

                       # Indexing and paths
                       self$index_of_non_terminal <- list()  # Index lookup for non-terminals
                       self$shortest_path <- list()  # Shortest paths in grammar
                     },

                     set_path = function(grammar_path) {
                       # Should validate input
                       if (!is.character(grammar_path) || length(grammar_path) != 1) {
                         stop("grammar_path must be a single character string")
                       }
                       self$grammar_file <- grammar_path
                       invisible(self)
                     },

                     set_pcfg_path = function(pcfg_path) {
                       # Should validate input
                       if (!is.character(pcfg_path) || length(pcfg_path) != 1) {
                         stop("pcfg_path must be a single character string")
                       }
                       self$pcfg_path <- pcfg_path
                       invisible(self)
                     },

                     set_min_init_tree_depth = function(min_tree_depth) {
                       self$max_init_depth <- min_tree_depth
                       invisible(self)
                     },

                     set_max_tree_depth = function(max_tree_depth) {
                       self$max_depth <- max_tree_depth
                       invisible(self)
                     },

                     get_max_depth = function() {
                       return(self$max_depth)
                     },

                     get_max_init_depth = function() {
                       return(self$max_init_depth)
                     },

                     read_grammar = function() {
                       if (is.null(self$grammar_file)) {
                         stop("You need to specify the path of the grammar file")
                       }

                       tryCatch({
                         # Read the grammar file line by line
                         lines <- readLines(self$grammar_file)

                         # Process each line
                         for (line in lines) {
                           # Skip comments and empty lines
                           if (!startsWith(line, "#") && nzchar(trimws(line))) {
                             if (grepl(self$PRODUCTION_SEPARATOR, line)) {
                               # Split into left side and productions
                               parts <- strsplit(line, self$RULE_SEPARATOR)[[1]]
                               if (length(parts) != 2) {
                                 stop(sprintf("Invalid grammar rule format in line: %s", line))
                               }

                               left_side <- trimws(parts[1])

                               # Validate left side is a non-terminal
                               if (!grepl(self$NT_PATTERN, left_side)) {
                                 stop("Left side not a non-terminal!")
                               }

                               # Add to non-terminals if not already present
                               if (!self$non_terminals$has(left_side)) {
                                 self$non_terminals$add(left_side)
                                 self$ordered_non_terminals$add(left_side)

                                 # Set start rule if not already set
                                 if (is.null(self$start_rule)) {
                                   self$start_rule <- list(left_side, self$NT)
                                 }
                               }

                               # Process productions
                               productions <- strsplit(parts[2], self$PRODUCTION_SEPARATOR)[[1]]
                               temp_productions <- list()

                               for (production in trimws(productions)) {
                                 temp_production <- list()

                                 if (!grepl(self$NT_PATTERN, production)) {
                                   # Handle terminal production
                                   if (production == "None") {
                                     production <- ""
                                   }
                                   if (!self$terminals$has(production)) {
                                     self$terminals$add(production)
                                   }
                                   temp_production[[length(temp_production) + 1]] <- list(production, self$T)
                                 } else {
                                   # Handle non-terminal production
                                   values <- unlist(regmatches(production,
                                                               gregexpr("<.+?>|[^<>]*", production)))

                                   for (value in values) {
                                     if (value != "") {
                                       if (!grepl(self$NT_PATTERN, value)) {
                                         sym <- list(value, self$T)
                                         if (!self$terminals$has(value)) {
                                           self$terminals$add(value)
                                         }
                                       } else {
                                         sym <- list(value, self$NT)
                                       }
                                       temp_production[[length(temp_production) + 1]] <- sym
                                     }
                                   }
                                 }
                                 temp_productions[[length(temp_productions) + 1]] <- temp_production
                               }

                               # Update max number of production rules
                               self$max_number_prod_rules <- max(self$max_number_prod_rules,
                                                                 length(temp_productions))

                               # Add to grammar if not already present
                               if (is.null(self$grammar[[left_side]])) {
                                 self$grammar[[left_side]] <- temp_productions
                               }
                             }
                           }
                         }

                         # Generate uniform PCFG
                         self$generate_uniform_pcfg()

                         # Load PCFG probabilities if path exists
                         if (!is.null(self$pcfg_path)) {
                           pcfg_data <- jsonlite::fromJSON(self$pcfg_path)
                           self$pcfg <- array(unlist(pcfg_data), dim = dim(pcfg_data))
                         }

                         # Find shortest path
                         self$find_shortest_path()

                       }, error = function(e) {
                         stop(sprintf("Error reading grammar file: %s", e$message))
                       })
                     },

                     find_shortest_path = function() {
                       # Initialize open_symbols as OrderedUniqueSet instead of list
                       open_symbols <- OrderedUniqueSet$new()

                       # Calculate minimum path for each non-terminal in grammar
                       for (nt in names(self$grammar)) {
                         depth <- self$minimum_path_calc(list(nt, 'NT'), open_symbols)
                       }
                     },

                     minimum_path_calc = function(current_symbol, open_symbols) {
                       # If current symbol is terminal, return 0
                       if (current_symbol[[2]] == self$T) {
                         return(0)
                       }

                       # Add current symbol value to open symbols using OrderedUniqueSet method
                       open_symbols$add(current_symbol[[1]])

                       # Process each derivation option for the current symbol
                       for (derivation_option in self$grammar[[current_symbol[[1]]]]) {
                         max_depth <- 0

                         # Initialize shortest path if not already present
                         key <- paste(current_symbol[[1]], current_symbol[[2]])
                         if (is.null(self$shortest_path[[key]])) {
                           self$shortest_path[[key]] <- list(999999)
                         }

                         # Check if any symbol in derivation_option is in open_symbols
                         has_open_symbol <- FALSE
                         for (sym in derivation_option) {
                           if (open_symbols$has(sym[[1]])) {  # Check symbol value only
                             has_open_symbol <- TRUE
                             break
                           }
                         }

                         # Skip this derivation if it contains open symbols
                         if (has_open_symbol) {
                           next
                         }

                         # Check if current_symbol appears in derivation_option
                         current_in_derivation <- FALSE
                         current_symbol_value <- current_symbol[[1]]  # Extract value once
                         for (sym in derivation_option) {
                           if (identical(sym[[1]], current_symbol_value)) {
                             current_in_derivation <- TRUE
                             break
                           }
                         }

                         # Process non-recursive derivations
                         if (!current_in_derivation) {
                           for (symbol in derivation_option) {
                             depth <- self$minimum_path_calc(symbol, open_symbols)
                             depth <- depth + 1
                             if (depth > max_depth) {
                               max_depth <- depth
                             }
                           }

                           # Update shortest path if new path is shorter or equal
                           if (max_depth < self$shortest_path[[key]][[1]]) {
                             # Found shorter path, reset and store new path
                             self$shortest_path[[key]] <- list(max_depth)
                             if (!derivation_option %in% self$shortest_path[[key]]) {
                               self$shortest_path[[key]][[length(self$shortest_path[[key]]) + 1]] <- derivation_option
                             }
                           } else if (max_depth == self$shortest_path[[key]][[1]]) {
                             # Found equal length path, add as alternative
                             if (!derivation_option %in% self$shortest_path[[key]]) {
                               self$shortest_path[[key]][[length(self$shortest_path[[key]]) + 1]] <- derivation_option
                             }
                           }
                         }
                       }

                       # Remove current symbol from open symbols using OrderedUniqueSet method
                       open_symbols$delete(current_symbol[[1]])  # Delete symbol value

                       # Return the shortest path length for this symbol
                       key <- paste(current_symbol[[1]], current_symbol[[2]])
                       return(self$shortest_path[[key]][[1]])
                     },

                     create_counter = function() {
                       # Initialize counter structure
                       self$counter <- list()

                       # For each key in the grammar, create a vector of zeros with length matching
                       # the number of production rules for that non-terminal
                       for (k in names(self$grammar)) {
                         self$counter[[k]] <- rep(0, length(self$grammar[[k]]))
                       }
                     },

                     generate_uniform_pcfg = function() {
                       if (length(self$grammar) == 0) {
                         stop("Cannot generate PCFG: grammar is empty")
                       }

                       # Create array with zeros
                       array <- matrix(0,
                                       nrow = length(self$grammar),
                                       ncol = self$max_number_prod_rules)

                       # For each non-terminal
                       for (i in seq_len(nrow(array))) {
                         nt <- names(self$grammar)[[i]]  # Use [[ for list element access
                         if (is.null(nt)) next  # Skip if name is NULL

                         number_probs <- length(self$grammar[[nt]])
                         if (number_probs > 0) {
                           # Calculate uniform probability
                           prob <- 1.0 / number_probs
                           # Assign probabilities to valid production rules
                           array[i, seq_len(number_probs)] <- prob
                         }

                         # Store index of non-terminal
                         if (is.null(self$index_of_non_terminal[[nt]])) {
                           self$index_of_non_terminal[[nt]] <- i
                         }
                       }

                       # Set PCFG and mask
                       self$pcfg <- array
                       self$pcfg_mask <- array != 0

                       invisible(self)  # For method chaining
                     },

                     generate_random_pcfg = function() {
                       if (length(self$grammar) == 0) {
                         stop("Cannot generate PCFG: grammar is empty")
                       }

                       # Create matrix with same dimensions as uniform PCFG
                       array <- matrix(0,
                                       nrow = length(self$grammar),
                                       ncol = self$max_number_prod_rules)

                       # For each non-terminal
                       for (i in seq_len(nrow(array))) {
                         nt <- names(self$grammar)[[i]]
                         if (is.null(nt)) next

                         number_probs <- length(self$grammar[[nt]])
                         if (number_probs > 0) {
                           # Generate random numbers from uniform distribution
                           random_values <- runif(number_probs)
                           # Normalize to ensure they sum to 1
                           normalized_probs <- random_values / sum(random_values)
                           # Assign to array
                           array[i, seq_len(number_probs)] <- normalized_probs
                         }

                         # Store index of non-terminal
                         if (is.null(self$index_of_non_terminal[[nt]])) {
                           self$index_of_non_terminal[[nt]] <- i
                         }
                       }

                       # Set PCFG and mask
                       self$pcfg <- array
                       self$pcfg_mask <- array != 0

                       invisible(self)
                     },

                     has_terminal = function(terminal) {
                       self$terminals$has(terminal)
                     },

                     has_non_terminal = function(non_terminal) {
                       self$non_terminals$has(non_terminal)
                     },

                     get_terminals = function() {
                       return(self$terminals)
                     },

                     get_non_terminals = function() {
                       return(self$ordered_non_terminals)
                     },

                     get_mask = function() {
                       if (is.null(self$pcfg_mask)) {
                         return(NULL)
                       }
                       return(self$pcfg_mask)
                     },

                     get_index_of_non_terminal = function() {
                       return(self$index_of_non_terminal)
                     },

                     count_number_of_options_in_production = function() {
                       if (is.null(self$number_of_options_by_non_terminal)) {
                         self$number_of_options_by_non_terminal <- list() # R list for Python dict
                         for (nt in self$get_non_terminals()$values()) { # Iterate using get_non_terminals() which now returns a list
                           self$number_of_options_by_non_terminal[[nt]] <- length(self$grammar[[nt]]) # length() in R for len() in Python
                         }
                       }
                       return(self$number_of_options_by_non_terminal)
                     },

                     compute_non_recursive_options = function() {
                       if (length(self$grammar) == 0) {
                         stop("Cannot compute options: grammar is empty")
                       }

                       # Initialize non_recursive_options if NULL
                       if (is.null(self$non_recursive_options)) {
                         self$non_recursive_options <- list()
                       }

                       # Iterate through each key in the grammar
                       for (key in names(self$grammar)) {
                         if (!self$non_terminals$has(key)) {
                           next
                         }

                         prob_non_recursive <- 0.0
                         non_recursive_prods <- list()

                         # Examine each production rule option
                         for (index in seq_along(self$grammar[[key]])) {
                           option <- self$grammar[[key]][[index]]

                           # Check if option is recursive
                           is_recursive <- FALSE
                           for (sym in option) {
                             if (!is.null(sym[[1]]) && identical(sym[[1]], key)) {
                               is_recursive <- TRUE
                               break
                             }
                           }

                           # If not recursive, process the option
                           if (!is_recursive) {
                             if (!is.null(self$index_of_non_terminal[[key]]) &&
                                 !is.null(self$pcfg[self$index_of_non_terminal[[key]], index])) {

                               prob_non_recursive <- prob_non_recursive +
                                 self$pcfg[self$index_of_non_terminal[[key]], index]

                               non_recursive_prods[[length(non_recursive_prods) + 1]] <-
                                 list(index = index, option = option)
                             }
                           }
                         }

                         # Store results
                         self$non_recursive_options[[key]] <- list(
                           productions = non_recursive_prods,
                           probability = round_prob(prob_non_recursive)
                         )
                       }

                       invisible(self)
                     },

                     list_non_recursive_productions = function(nt) {
                       if (!self$non_terminals$has(nt)) {
                         stop("Non-terminal not found in grammar")
                       }

                       non_recursive_elements <- list()
                       for (options in self$grammar[[nt]]) {
                         is_recursive <- FALSE
                         for (option in options) {
                           if (!self$terminals$has(option[[1]]) && identical(option[[1]], nt)) {
                             is_recursive <- TRUE
                             break
                           }
                         }
                         if (!is_recursive) {
                           non_recursive_elements[[length(non_recursive_elements) + 1]] <- options
                         }
                       }
                       return(non_recursive_elements)
                     },

                     recursive_individual_creation = function(genome, symbol, current_depth) {
                       # Generate random uniform value between 0 and 1
                       codon <- runif(1)

                       # Get index of the non-terminal in PCFG
                       if (!self$non_terminals$has(symbol)) {
                         stop(sprintf("Symbol '%s' is not a non-terminal in the grammar", symbol))
                       }
                       nt_index <- self$index_of_non_terminal[[symbol]]

                       # Handle depth constraints and choose expansion possibility
                       if (current_depth > self$max_init_depth) {
                         # Get shortest path for the current symbol
                         shortest_path_key <- paste(symbol, 'NT')
                         if (is.null(self$shortest_path[[shortest_path_key]])) {
                           stop(sprintf("No shortest path found for symbol '%s'", symbol))
                         }
                         shortest_path <- self$shortest_path[[shortest_path_key]]

                         # Calculate probability for non-recursive options
                         prob_non_recursive <- 0.0
                         for (rule in shortest_path[-1]) {  # Skip first element which is the depth
                           index <- which(sapply(self$grammar[[symbol]], function(x) identical(x, rule)))
                           if (length(index) > 0) {
                             prob_non_recursive <- prob_non_recursive + self$pcfg[nt_index, index]
                           }
                         }

                         # Choose expansion based on normalized probabilities
                         prob_aux <- 0.0
                         expansion_possibility <- NULL
                         for (rule in shortest_path[-1]) {
                           index <- which(sapply(self$grammar[[symbol]], function(x) identical(x, rule)))
                           if (length(index) > 0) {
                             new_prob <- self$pcfg[nt_index, index] / prob_non_recursive
                             prob_aux <- prob_aux + new_prob
                             if (codon <= round(prob_aux, 3)) {
                               expansion_possibility <- index
                               break
                             }
                           }
                         }
                       } else {
                         # Standard probability-based selection
                         prob_aux <- 0.0
                         expansion_possibility <- NULL
                         for (index in seq_along(self$grammar[[symbol]])) {
                           prob_aux <- prob_aux + self$pcfg[nt_index, index]
                           if (codon <= round(prob_aux, 3)) {
                             expansion_possibility <- index
                             break
                           }
                         }
                       }

                       # Add chosen expansion to genome
                       non_terminals <- self$get_non_terminals()$values()
                       nt_position <- which(sapply(non_terminals, function(x) identical(x, symbol)))
                       if (length(nt_position) == 0) {
                         stop(sprintf("Position not found for symbol '%s' in non-terminals", symbol))
                       }
                       genome[[nt_position]] <- append(genome[[nt_position]],
                                                       list(c(expansion_possibility, codon, current_depth)))

                       # Get expansion symbols
                       expansion_symbols <- self$grammar[[symbol]][[expansion_possibility]]

                       # Track depths for recursive calls
                       depths <- c(current_depth)

                       # Recursively process non-terminal symbols
                       for (sym in expansion_symbols) {
                         if (!self$terminals$has(sym[[1]])) {  # Use has() to check if it's not a terminal
                           depths <- c(depths,
                                       self$recursive_individual_creation(genome, sym[[1]], current_depth + 1))
                         }
                       }

                       # Return maximum depth reached
                       return(max(depths))
                     },

                     get_probabilities_non_terminal = function(grammar, nt_index) {
                       # Input validation
                       if (is.null(grammar) || is.null(nt_index)) {
                         return(NULL)
                       }
                       if (!is.numeric(nt_index)) {
                         stop("nt_index must be numeric")
                       }
                       if (!(is.matrix(grammar) || is.array(grammar))) {
                         stop("grammar must be a matrix or array")
                       }

                       # Return the entire row of probabilities
                       tryCatch({
                         if (is.matrix(grammar)) {
                           return(grammar[nt_index,])
                         } else if (is.array(grammar)) {
                           return(grammar[nt_index,])
                         }
                       }, error = function(e) {
                         stop(sprintf("Error accessing probabilities: %s", e$message))
                       })
                     },

                     get_probability = function(grammar, nt_index, index) {
                       # Input validation
                       if (is.null(grammar) || is.null(nt_index) || is.null(index)) {
                         return(NULL)
                       }
                       if (!is.numeric(nt_index) || !is.numeric(index)) {
                         stop("nt_index and index must be numeric")
                       }
                       if (!(is.matrix(grammar) || is.array(grammar))) {
                         stop("grammar must be a matrix or array")
                       }

                       # Handle both matrix and array cases
                       tryCatch({
                         if (is.matrix(grammar)) {
                           return(grammar[nt_index, index])
                         } else if (is.array(grammar)) {
                           return(grammar[nt_index, index])
                         }
                       }, error = function(e) {
                         stop(sprintf("Error accessing probability: %s", e$message))
                       })
                     },

                     mapping = function(mapping_rules, positions_to_map = NULL, needs_r_filter = FALSE) {
                       # Initialize positions_to_map if NULL
                       if (is.null(positions_to_map)) {
                         positions_to_map <- rep(0, length(self$ordered_non_terminals$values()))
                       }

                       # Initialize output vector to store text parts
                       output <- character(0)

                       # Call recursive_mapping to generate the output
                       max_depth <- self$recursive_mapping(
                         mapping_rules = mapping_rules,
                         positions_to_map = positions_to_map,
                         current_sym = self$start_rule,
                         current_depth = 0,
                         output = output
                       )

                       # Combine output into single string
                       output <- paste(output, collapse = "")

                       # Apply R filter if needed
                       if (endsWith(self$grammar_file, "pybnf")) {
                         output <- r_filter(output, needs_r_filter)
                       }

                       # Return both output and max_depth as a list
                       return(list(
                         output = output,
                         max_depth = max_depth
                       ))
                     },

                     recursive_mapping = function(mapping_rules, positions_to_map, current_sym, current_depth, output) {
                       # Validate inputs
                       if (!is.list(mapping_rules) || !is.numeric(current_depth) || !is.character(output)) {
                         stop("Invalid input types for recursive_mapping")
                       }

                       # Initialize depths vector with current depth
                       depths <- c(current_depth)

                       if (self$terminals$has(current_sym[[1]])) {
                         # If current symbol is terminal, append to output
                         output[length(output) + 1] <- current_sym[[1]]
                       } else {
                         # Validate current symbol is a non-terminal
                         if (!self$non_terminals$has(current_sym[[1]])) {
                           stop(sprintf("Unknown non-terminal symbol: %s", current_sym[[1]]))
                         }

                         # Get position of current non-terminal in ordered list
                         current_sym_pos <- which(sapply(self$ordered_non_terminals$values(),
                                                         function(x) identical(x, current_sym[[1]])))
                         if (length(current_sym_pos) == 0) {
                           stop(sprintf("Symbol %s not found in ordered non-terminals", current_sym[[1]]))
                         }

                         # Get possible choices for this non-terminal
                         choices <- self$grammar[[current_sym[[1]]]]
                         if (is.null(choices)) {
                           stop(sprintf("No production rules found for symbol: %s", current_sym[[1]]))
                         }

                         # Get shortest path info
                         shortest_path_key <- paste(current_sym[[1]], current_sym[[2]])
                         shortest_path <- self$shortest_path[[shortest_path_key]]
                         if (is.null(shortest_path)) {
                           stop(sprintf("No shortest path found for symbol: %s", current_sym[[1]]))
                         }

                         # Get index in PCFG
                         nt_index <- self$index_of_non_terminal[[current_sym[[1]]]]
                         if (is.null(nt_index)) {
                           stop(sprintf("No PCFG index found for symbol: %s", current_sym[[1]]))
                         }

                         # Handle expansion possibility selection
                         if (positions_to_map[current_sym_pos] >= length(mapping_rules[[current_sym_pos]])) {
                           # Generate new expansion if we've exhausted mapping rules
                           codon <- runif(1)

                           if (current_depth > self$max_depth) {
                             # Calculate probability for non-recursive options
                             prob_non_recursive <- 0.0
                             for (rule in shortest_path[-1]) {
                               index <- which(sapply(self$grammar[[current_sym[[1]]]],
                                                     function(x) identical(x, rule)))
                               if (length(index) > 0) {
                                 prob_non_recursive <- prob_non_recursive + self$pcfg[nt_index, index]
                               }
                             }

                             # Choose expansion based on normalized probabilities
                             prob_aux <- 0.0
                             expansion_possibility <- NULL
                             for (rule in shortest_path[-1]) {
                               index <- which(sapply(self$grammar[[current_sym[[1]]]],
                                                     function(x) identical(x, rule)))
                               if (length(index) > 0) {
                                 new_prob <- self$pcfg[nt_index, index] / prob_non_recursive
                                 prob_aux <- prob_aux + new_prob
                                 if (codon <= round(prob_aux, 3)) {
                                   expansion_possibility <- index
                                   break
                                 }
                               }
                             }
                           } else {
                             # Standard probability-based selection
                             prob_aux <- 0.0
                             expansion_possibility <- NULL
                             for (index in seq_along(self$grammar[[current_sym[[1]]]])) {
                               prob_aux <- prob_aux + self$pcfg[nt_index, index]
                               if (codon <= round(prob_aux, 3)) {
                                 expansion_possibility <- index
                                 break
                               }
                             }
                           }

                           # Add new mapping rule
                           mapping_rules[[current_sym_pos]] <- append(
                             mapping_rules[[current_sym_pos]],
                             list(c(expansion_possibility, codon, current_depth))
                           )
                         } else {
                           # Re-use existing mapping rule with updated probabilities
                           codon <- mapping_rules[[current_sym_pos]][[positions_to_map[current_sym_pos] + 1]][2]

                           if (current_depth > self$max_depth) {
                             # Calculate probability for non-recursive options
                             prob_non_recursive <- 0.0
                             for (rule in shortest_path[-1]) {
                               index <- which(sapply(self$grammar[[current_sym[[1]]]],
                                                     function(x) identical(x, rule)))
                               if (length(index) > 0) {
                                 prob_non_recursive <- prob_non_recursive + self$pcfg[nt_index, index]
                               }
                             }

                             # Choose expansion based on normalized probabilities
                             prob_aux <- 0.0
                             for (rule in shortest_path[-1]) {
                               index <- which(sapply(self$grammar[[current_sym[[1]]]],
                                                     function(x) identical(x, rule)))
                               if (length(index) > 0) {
                                 new_prob <- self$pcfg[nt_index, index] / prob_non_recursive
                                 prob_aux <- prob_aux + new_prob
                                 if (codon <= round(prob_aux, 3)) {
                                   expansion_possibility <- index
                                   break
                                 }
                               }
                             }
                           } else {
                             # Standard probability-based selection
                             prob_aux <- 0.0
                             for (index in seq_along(self$grammar[[current_sym[[1]]]])) {
                               prob_aux <- prob_aux + self$pcfg[nt_index, index]
                               if (codon <= round(prob_aux, 3)) {
                                 expansion_possibility <- index
                                 break
                               }
                             }
                           }
                         }

                         # Update mapping rules with new expansion
                         mapping_rules[[current_sym_pos]][[positions_to_map[current_sym_pos] + 1]] <-
                           c(expansion_possibility, codon, current_depth)

                         # Get current production and increment position
                         current_production <- expansion_possibility
                         positions_to_map[current_sym_pos] <- positions_to_map[current_sym_pos] + 1

                         # Get next symbols to expand
                         next_to_expand <- choices[[current_production]]

                         # Recursively process each symbol
                         for (next_sym in next_to_expand) {
                           depths <- c(depths, self$recursive_mapping(
                             mapping_rules = mapping_rules,
                             positions_to_map = positions_to_map,
                             current_sym = next_sym,
                             current_depth = current_depth + 1,
                             output = output
                           ))
                         }
                       }

                       # Return maximum depth reached
                       return(max(depths))
                     },

                     get_non_recursive_options = function(symbol) {
                       if (is.null(symbol)) {
                         stop("Symbol cannot be NULL")
                       }
                       if (!self$non_terminals$has(symbol)) {
                         stop(sprintf("Symbol '%s' not found in non-terminals", symbol))
                       }
                       return(self$non_recursive_options[[symbol]])
                     },

                     get_dict = function() {
                       # def get_dict(self):
                       return(self$grammar)
                     },

                     get_pcfg = function() {
                       # def get_pcfg(self):
                       return(self$pcfg)
                     },

                     get_shortest_path = function() {
                       # def get_shortest_path(self):
                       return(self$shortest_path)
                     },

                     get_start_rule = function() {
                       # def get_start_rule(self):
                       return(self$start_rule)
                     },

                     find_common_terminals = function(other_grammar) {
                       # Validate input
                       if (!inherits(other_grammar, "Grammar")) {
                         stop("other_grammar must be a Grammar object")
                       }

                       # Get terminals from both grammars using the values() method
                       my_terminals <- self$terminals$values()
                       other_terminals <- other_grammar$terminals$values()

                       # Create new OrderedUniqueSet for common terminals
                       common <- OrderedUniqueSet$new()

                       # Add terminals that exist in both grammars
                       for (term in my_terminals) {
                         if (other_grammar$terminals$has(term)) {
                           common$add(term)
                         }
                       }

                       return(common)
                     },

                     vector_to_ordered_set = function(vec) {
                       set <- OrderedUniqueSet$new()
                       for (elem in vec) {
                         set$add(elem)
                       }
                       return(set)
                     },

                     ordered_set_to_vector = function(set) {
                       return(unlist(set$values()))
                     },

                     print = function() {
                       # String representation of the grammar
                       grammar <- self$grammar
                       text <- ""

                       # Get non-terminals from OrderedUniqueSet
                       non_terminals <- self$get_non_terminals()$values()

                       # Process each non-terminal
                       for (key in non_terminals) {
                         # Add non-terminal and separator
                         text <- paste0(text, key, " ::= ")

                         # Get production rules for this non-terminal
                         options_list <- grammar[[key]]

                         # Process each production rule
                         for (i in seq_along(options_list)) {
                           options <- options_list[[i]]

                           # Process each symbol in the production
                           for (option in options) {
                             text <- paste0(text, option[[1]])
                           }

                           # Add separator if not the last production
                           if (i < length(options_list)) {
                             text <- paste0(text, " | ")
                           }
                         }

                         # Add newline after each non-terminal's rules
                         text <- paste0(text, "\n")
                       }

                       # Print the grammar representation
                       cat(text)

                       # Return invisibly for method chaining
                       invisible(self)
                     }
                   )
)

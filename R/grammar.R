# Helper functions
round_prob <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0) {
    return(x)
  }
  if (any(is.na(x))) {
    return(x)
  } # Preserve NA values
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

    tryCatch(
      {
        # Attempt to style the text as R code using styler::style_text
        formatted_txt <- styler::style_text(txt_no_markers)
        # styler::style_text returns a list of formatted code lines.
        txt <- paste(formatted_txt, collapse = "\n") # Convert to single string
      },
      error = function(e) {
        # If styler fails to parse or style (e.g., if 'txt' is not valid R-like code)
        warning("styler formatting failed: ", e$message, ". Returning unstyled text (markers removed).")
        # In case of error, return the original text without styling
        return(txt_no_markers)
      }
    )
  }
  return(trimws(txt))
}

Grammar <- R6::R6Class(
  "Grammar",
  public = list(
    NT = "NT",
    T = "T",
    NT_PATTERN = "(<[^>]+>)",
    RULE_SEPARATOR = "::=",
    PRODUCTION_SEPARATOR = "|",
    grammar = list(),
    non_terminals = NULL,
    terminals = NULL,
    ordered_non_terminals = NULL,
    pcfg = NULL,
    pcfg_mask = NULL,
    max_depth = NULL,
    max_init_depth = NULL,
    grammar_file = NULL,
    start_rule = NULL,
    index_of_non_terminal = list(),
    non_recursive_options = list(),
    shortest_path = list(),
    initialize = function() {
      self$non_terminals <- OrderedUniqueSet$new()
      self$terminals <- OrderedUniqueSet$new()
      self$ordered_non_terminals <- OrderedUniqueSet$new()
    },
    set_path = function(grammar_path) {
      self$grammar_file <- grammar_path
      invisible(self)
    },
    set_max_tree_depth = function(max_tree_depth) {
      self$max_depth <- max_tree_depth
      invisible(self)
    },
    set_min_init_tree_depth = function(min_tree_depth) {
      self$max_init_depth <- min_tree_depth
      invisible(self)
    },
    parse_production = function(production) {
      tokens <- stringr::str_match_all(production, "<[^>]+>|[^<>]+")[[1]][,1]
      tokens <- trimws(tokens)

      result <- list()
      for (token in tokens) {
        if (grepl(self$NT_PATTERN, token)) {
          result[[length(result) + 1]] <- list(token, self$NT)
          if (!self$non_terminals$has(token)) {
            self$non_terminals$add(token)
            self$ordered_non_terminals$add(token)
          }
        } else if (nchar(token) > 0) {
          result[[length(result) + 1]] <- list(token, self$T)
          if (!self$terminals$has(token)) {
            self$terminals$add(token)
          }
        }
      }
      result
    },
    read_grammar = function() {
      if (is.null(self$grammar_file)) {
        stop("Grammar file path not set")
      }

      lines <- readLines(self$grammar_file)
      self$grammar <- list()

      for (line in lines) {
        line <- trimws(line)
        if (nchar(line) == 0 || startsWith(line, "#")) next

        parts <- strsplit(line, self$RULE_SEPARATOR)[[1]]
        if (length(parts) != 2) next

        lhs <- trimws(parts[1])
        rhs <- trimws(parts[2])

        productions <- strsplit(rhs, paste0("\\s*\\", self$PRODUCTION_SEPARATOR, "\\s*"))[[1]]
        temp_productions <- lapply(productions, self$parse_production)

        self$grammar[[lhs]] <- temp_productions
      }

      if (length(self$grammar) > 0 && is.null(self$start_rule)) {
        first_nt <- names(self$grammar)[1]
        self$start_rule <- list(first_nt, self$NT)
      }

      # Initialize index lookup for non-terminals
      self$index_of_non_terminal <- list()
      for (i in seq_along(self$grammar)) {
        nt <- names(self$grammar)[i]
        self$index_of_non_terminal[[nt]] <- i
      }

      # Initialize supporting structures
      self$generate_uniform_pcfg()
      self$compute_non_recursive_options()
      self$find_shortest_path()

      invisible(self)
    },
    generate_uniform_pcfg = function() {
      if (length(self$grammar) == 0) return()

      max_rules <- max(sapply(self$grammar, length))
      self$pcfg <- matrix(0,
                          nrow = length(self$grammar),
                          ncol = max_rules)
      rownames(self$pcfg) <- names(self$grammar)

      for (i in seq_along(self$grammar)) {
        nt <- names(self$grammar)[i]
        n_rules <- length(self$grammar[[nt]])
        if (n_rules > 0) {
          self$pcfg[i, 1:n_rules] <- 1/n_rules
        }
      }

      self$pcfg_mask <- self$pcfg != 0
    },
    update_pcfg = function(best_individual, learning_factor) {
      if (is.null(best_individual$genotype)) {
        stop("Invalid individual: no genotype found")
      }

      for (nt in names(best_individual$genotype)) {
        nt_index <- self$index_of_non_terminal[[nt]]
        if (is.null(nt_index)) next

        # Count rule usage
        rule_counts <- integer(ncol(self$pcfg))
        rule_usage <- best_individual$genotype[[nt]]

        for (gene in rule_usage) {
          if (!is.null(gene) && length(gene) > 0) {
            rule_counts[gene[1]] <- rule_counts[gene[1]] + 1
          }
        }

        total_used <- sum(rule_counts)
        if (total_used == 0) next

        current_probs <- self$pcfg[nt_index,]
        new_probs <- current_probs

        # Update probabilities using equations from paper
        for (j in seq_along(rule_counts)) {
          if (rule_counts[j] > 0) {
            # Increase probability for used rules
            new_probs[j] <- min(
              current_probs[j] + learning_factor * (rule_counts[j] / total_used),
              1.0
            )
          } else {
            # Decrease probability for unused rules
            new_probs[j] <- max(
              current_probs[j] - learning_factor * current_probs[j],
              0.0
            )
          }
        }

        # Normalize the row
        self$pcfg[nt_index,] <- private$normalize_row_probabilities(new_probs)
      }
    },
    verify_pcfg = function() {
      if (is.null(self$pcfg)) return(FALSE)

      for (i in seq_len(nrow(self$pcfg))) {
        active_probs <- self$pcfg[i, self$pcfg_mask[i,]]
        if (length(active_probs) > 0) {
          row_sum <- sum(active_probs)
          if (abs(row_sum - 1) > 1e-6) {
            warning(sprintf(
              "Row %d probabilities sum to %f, normalizing",
              i, row_sum
            ))
            self$pcfg[i, self$pcfg_mask[i,]] <- active_probs / row_sum
          }
        }
      }
      TRUE
    },
    find_shortest_path = function() {
      open_symbols <- OrderedUniqueSet$new()

      for (nt in names(self$grammar)) {
        depth <- self$minimum_path_calc(list(nt, self$NT), open_symbols)
      }

      invisible(self)
    },
    minimum_path_calc = function(current_symbol, open_symbols) {
      if (current_symbol[[2]] == self$T) {
        return(0)
      }

      open_symbols$add(current_symbol[[1]])

      key <- paste(current_symbol[[1]], current_symbol[[2]])
      if (is.null(self$shortest_path[[key]])) {
        self$shortest_path[[key]] <- list(depth = 999999, options = list())
      }

      for (derivation_option in self$grammar[[current_symbol[[1]]]]) {
        max_depth <- 0

        # Check for open symbols
        has_open_symbol <- FALSE
        for (sym in derivation_option) {
          if (open_symbols$has(sym[[1]])) {
            has_open_symbol <- TRUE
            break
          }
        }
        if (has_open_symbol) next

        # Check for recursion
        current_in_derivation <- FALSE
        current_symbol_value <- current_symbol[[1]]
        for (sym in derivation_option) {
          if (identical(sym[[1]], current_symbol_value)) {
            current_in_derivation <- TRUE
            break
          }
        }

        if (!current_in_derivation) {
          for (symbol in derivation_option) {
            depth <- self$minimum_path_calc(symbol, open_symbols)
            depth <- depth + 1
            if (depth > max_depth) {
              max_depth <- depth
            }
          }

          # Compare with current shortest path
          if (max_depth < self$shortest_path[[key]]$depth) {
            self$shortest_path[[key]]$depth <- max_depth
            self$shortest_path[[key]]$options <- list(derivation_option)
          } else if (max_depth == self$shortest_path[[key]]$depth) {
            # Check if this option is already present
            is_duplicate <- FALSE
            for (existing_option in self$shortest_path[[key]]$options) {
              if (identical(existing_option, derivation_option)) {
                is_duplicate <- TRUE
                break
              }
            }
            if (!is_duplicate) {
              self$shortest_path[[key]]$options <- c(
                self$shortest_path[[key]]$options,
                list(derivation_option)
              )
            }
          }
        }
      }

      open_symbols$delete(current_symbol[[1]])
      return(self$shortest_path[[key]]$depth)
    },
    compute_non_recursive_options = function() {
      if (length(self$grammar) == 0) {
        stop("Cannot compute options: grammar is empty")
      }

      for (nt in names(self$grammar)) {
        if (!self$non_terminals$has(nt)) next

        prob_non_recursive <- 0.0
        non_recursive_prods <- list()

        for (index in seq_along(self$grammar[[nt]])) {
          option <- self$grammar[[nt]][[index]]

          is_recursive <- FALSE
          for (sym in option) {
            if (!is.null(sym[[1]]) && identical(sym[[1]], nt)) {
              is_recursive <- TRUE
              break
            }
          }

          if (!is_recursive) {
            nt_index <- self$index_of_non_terminal[[nt]]
            if (!is.null(nt_index) && !is.null(self$pcfg[nt_index, index])) {
              prob_non_recursive <- prob_non_recursive + self$pcfg[nt_index, index]
              non_recursive_prods[[length(non_recursive_prods) + 1]] <- list(
                index = index,
                option = option
              )
            }
          }
        }

        self$non_recursive_options[[nt]] <- list(
          productions = non_recursive_prods,
          probability = round(prob_non_recursive, 3)
        )
      }

      invisible(self)
    },
    recursive_individual_creation = function(genotype, symbol, current_depth = 0) {
      # Initialize if empty
      if (length(genotype) == 0) {
        for (nt in self$ordered_non_terminals$values()) {
          genotype[[nt[[1]]]] <- list()
        }
      }

      codon <- runif(1)
      nt_index <- self$index_of_non_terminal[[symbol]]

      # Select rule using PCFG probabilities
      prob_aux <- 0
      expansion_possibility <- 1
      for(index in seq_along(self$grammar[[symbol]])) {
        prob_aux <- prob_aux + self$pcfg[nt_index, index]
        if(codon <= round(prob_aux, 3)) {
          expansion_possibility <- index
          break
        }
      }

      # Add this choice to the genotype
      genotype[[symbol]][[length(genotype[[symbol]]) + 1]] <-
        c(expansion_possibility, codon, current_depth)

      # Recursively handle expansion
      expansion_symbols <- self$grammar[[symbol]][[expansion_possibility]]
      depths <- current_depth

      for(sym in expansion_symbols) {
        if(sym[[2]] == self$NT) {
          result <- self$recursive_individual_creation(
            genotype,
            sym[[1]],
            current_depth + 1
          )
          genotype <- result$genotype
          depths <- c(depths, result$depth)
        }
      }

      list(
        genotype = genotype,
        depth = max(depths)
      )
    },
    mapping = function(mapping_rules, positions_to_map = NULL) {
      if (is.null(positions_to_map)) {
        positions_to_map <- rep(0, length(self$ordered_non_terminals$values()))
      }

      result <- self$recursive_mapping(
        mapping_rules,
        positions_to_map,
        self$start_rule,
        0,
        character(0)
      )

      list(
        phenotype = paste(result$output, collapse = ""),
        mapping_result = result$depth,
        positions = result$positions
      )
    },
    recursive_mapping = function(mapping_rules, positions_to_map, current_sym, current_depth, output) {
      depths <- current_depth

      if (identical(current_sym[[2]], self$T)) {
        output <- c(output, current_sym[[1]])
      } else {
        current_sym_pos <- which(sapply(self$ordered_non_terminals$values(),
                                        function(x) identical(x, current_sym[[1]])))

        rule_pos <- positions_to_map[current_sym_pos] + 1
        expansion_possibility <- mapping_rules[[current_sym_pos]][[rule_pos]][[1]]
        positions_to_map[current_sym_pos] <- rule_pos

        next_symbols <- self$grammar[[current_sym[[1]]]][[expansion_possibility]]

        for(next_sym in next_symbols) {
          result <- self$recursive_mapping(
            mapping_rules,
            positions_to_map,
            next_sym,
            current_depth + 1,
            output
          )
          output <- result$output
          depths <- c(depths, result$depth)
          positions_to_map <- result$positions
        }
      }

      list(
        output = output,
        depth = max(depths),
        positions = positions_to_map
      )
    },
    has_terminal = function(terminal) {
      self$terminals$has(terminal)
    },
    has_non_terminal = function(non_terminal) {
      self$non_terminals$has(non_terminal)
    },
    get_terminals = function() {
      self$terminals
    },
    get_non_terminals = function() {
      self$ordered_non_terminals
    },
    get_pcfg = function() {
      if (!self$verify_pcfg()) {
        stop("PCFG verification failed")
      }
      self$pcfg
    },
    get_mask = function() {
      if (is.null(self$pcfg)) {
        return(NULL)
      }
      self$pcfg != 0
    },
    get_index_of_non_terminal = function() {
      self$index_of_non_terminal
    },
    get_non_recursive_options = function(symbol) {
      if (is.null(symbol)) stop("Symbol cannot be NULL")
      if (!self$non_terminals$has(symbol)) {
        stop(sprintf("Symbol '%s' not found in non-terminals", symbol))
      }
      self$non_recursive_options[[symbol]]
    },
    get_dict = function() {
      self$grammar
    },
    get_shortest_path = function() {
      self$shortest_path
    },
    get_start_rule = function() {
      self$start_rule
    },
    get_max_depth = function() {
      self$max_depth
    },
    get_max_init_depth = function() {
      self$max_init_depth
    },
    normalize_row_probabilities = function(probs) {
      # Handle zero-sum case
      if (sum(probs) == 0) {
        n <- length(probs)
        return(rep(1/n, n))
      }

      # Normal normalization
      probs / sum(probs)
    }
  )
)

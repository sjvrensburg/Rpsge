#' Grammar Analysis for Diverse Initialization (GADI Phase 1)
#'
#' Analyzes a grammar to identify key characteristics that influence diversity
#' potential, including recursive structures, rule diversity indices, and non-terminal
#' dependencies. This analysis forms the foundation of the Grammar-Aware Diverse
#' Initialization (GADI) method.
#'
#' @param grammar Grammar structure from read_pcfg_grammar()
#' @param use_semantic_awareness Boolean indicating whether to use semantic-aware
#'   diversity metrics (requires a semantic grouping function)
#' @param semantic_group_fn Optional function that groups semantically equivalent
#'   rules. This function should take two arguments:
#'   \itemize{
#'     \item \code{rules}: A list of production rules for a specific non-terminal
#'     \item \code{non_terminal}: The name of the non-terminal these rules belong to
#'   }
#'   The function should return a list of rule groups, where each group is a list
#'   containing rules that are considered semantically equivalent. For example, if
#'   a non-terminal has 5 rules and rules 1 and 3 are equivalent, while rules 2, 4,
#'   and 5 are each distinct, the function should return a list with 4 elements:
#'   \code{list(list(rule1, rule3), list(rule2), list(rule4), list(rule5))}.
#'   This parameter is required if use_semantic_awareness is TRUE.
#'
#' @return A list containing grammar analysis information:
#'   \itemize{
#'     \item \code{start_symbol}: The grammar's start symbol
#'     \item \code{recursive_nts}: Non-terminals with recursive rules
#'     \item \code{dependencies}: List of dependencies between non-terminals
#'     \item \code{rule_counts}: Number of rules for each non-terminal
#'     \item \code{non_terminal_depth}: Minimum tree depths for each non-terminal
#'     \item \code{rule_diversity_index}: Rule Diversity Index (RDI) values
#'     \item \code{semantic_rule_diversity_index}: Semantic-aware RDI values
#'     \item \code{sorted_by_diversity}: Non-terminals sorted by diversity potential
#'     \item \code{sorted_by_rules}: Non-terminals sorted by rule count
#'     \item \code{high_diversity_paths}: Promising grammar paths for diversity
#'   }
#'
#' @examples
#' # Simple grammar analysis without semantic awareness
#' grammar <- read_pcfg_grammar(system.file("extdata", "expr.bnf", package = "Rpsge"))
#' analysis <- analyze_grammar_for_diversity(grammar)
#'
#' # Grammar analysis with semantic awareness for mathematical expressions
#' # Create a semantic grouping function for commutative operations
#' semantic_grouper <- function(rules, non_terminal) {
#'   # Initialize groups
#'   groups <- list()
#'
#'   # Track which rules have been grouped
#'   grouped <- logical(length(rules))
#'
#'   # For each rule
#'   for (i in seq_along(rules)) {
#'     if (grouped[i]) next
#'
#'     current_group <- list(rules[[i]])
#'     grouped[i] <- TRUE
#'
#'     # Check if this is a binary operation with a commutative operator
#'     if (length(rules[[i]]$symbols) == 3 &&
#'         rules[[i]]$symbols[2] %in% c("+", "*")) {
#'
#'       # Look for rules with the same operator but operands swapped
#'       for (j in seq_along(rules)) {
#'         if (grouped[j] || i == j) next
#'
#'         # Check for commutativity pattern
#'         if (length(rules[[j]]$symbols) == 3 &&
#'             rules[[i]]$symbols[2] == rules[[j]]$symbols[2] &&
#'             rules[[i]]$symbols[1] == rules[[j]]$symbols[3] &&
#'             rules[[i]]$symbols[3] == rules[[j]]$symbols[1] &&
#'             all(rules[[i]]$types == rules[[j]]$types)) {
#'
#'           current_group <- c(current_group, list(rules[[j]]))
#'           grouped[j] <- TRUE
#'         }
#'       }
#'     }
#'
#'     groups <- c(groups, list(current_group))
#'   }
#'
#'   return(groups)
#' }
#'
#' # Analyze with semantic awareness
#' semantic_analysis <- analyze_grammar_for_diversity(
#'   grammar,
#'   use_semantic_awareness = TRUE,
#'   semantic_group_fn = semantic_grouper
#' )
#'
#' @export
analyze_grammar_for_diversity <- function(grammar,
                                          use_semantic_awareness = FALSE,
                                          semantic_group_fn = NULL) {
  # Validate inputs
  if (!is.list(grammar) || !all(c("non_terminals", "terminals", "rules", "start_symbol") %in% names(grammar))) {
    stop("Invalid grammar structure. Must contain non_terminals, terminals, rules, and start_symbol.")
  }

  if (use_semantic_awareness && is.null(semantic_group_fn)) {
    stop("Semantic-aware analysis requires a semantic_group_fn function")
  }

  # Initialize analysis results
  analysis <- list(
    start_symbol = grammar$start_symbol,
    recursive_nts = character(0),
    dependencies = list(),
    rule_counts = list(),
    non_terminal_depth = list(),
    rule_diversity_index = list(),
    semantic_rule_diversity_index = list(),
    sorted_by_diversity = character(0),
    sorted_by_rules = character(0),
    high_diversity_paths = list()
  )

  # 1. Identify recursive rules using the package's existing function
  grammar_with_flags <- identify_recursive_rules(grammar)

  # Extract recursive non-terminals from the grammar with flags
  for (nt in names(grammar_with_flags$rules)) {
    for (i in seq_along(grammar_with_flags$rules[[nt]])) {
      if (!is.null(grammar_with_flags$rules[[nt]][[i]]$is_recursive) &&
        grammar_with_flags$rules[[nt]][[i]]$is_recursive) {
        if (!(nt %in% analysis$recursive_nts)) {
          analysis$recursive_nts <- c(analysis$recursive_nts, nt)
        }
      }
    }
  }

  # 2. Identify dependencies between non-terminals
  for (nt in grammar$non_terminals) {
    # Initialize dependencies for this non-terminal
    analysis$dependencies[[nt]] <- character(0)

    # Count rules for this non-terminal
    if (nt %in% names(grammar$rules)) {
      analysis$rule_counts[[nt]] <- length(grammar$rules[[nt]])

      # Check each rule for dependencies
      for (rule_idx in seq_along(grammar$rules[[nt]])) {
        rule <- grammar$rules[[nt]][[rule_idx]]

        # Look for non-terminal symbols in this rule
        nt_indices <- which(rule$types == "NT")
        if (length(nt_indices) > 0) {
          dependent_nts <- unique(rule$symbols[nt_indices])

          # Add to dependencies if not already included
          for (dep_nt in dependent_nts) {
            if (!(dep_nt %in% analysis$dependencies[[nt]])) {
              analysis$dependencies[[nt]] <- c(analysis$dependencies[[nt]], dep_nt)
            }
          }
        }
      }
    } else {
      analysis$rule_counts[[nt]] <- 0
    }
  }

  # 3. Calculate minimum tree depths using the package's existing function
  analysis$non_terminal_depth <- calculate_min_depths(grammar)

  # 4. Calculate Rule Diversity Index (RDI) for each non-terminal
  for (nt in grammar$non_terminals) {
    # Skip if non-terminal has no rules defined
    if (!(nt %in% names(grammar$rules))) {
      analysis$rule_diversity_index[[nt]] <- 0
      next
    }

    # Get the number of rules for this non-terminal
    rule_count <- length(grammar$rules[[nt]])

    if (rule_count <= 1) {
      # Non-terminals with 0 or 1 rules have 0 diversity
      analysis$rule_diversity_index[[nt]] <- 0
      next
    }

    # Calculate RDI as rule_count * log(rule_count)
    # This follows the formula defined in the GADI paper
    analysis$rule_diversity_index[[nt]] <- rule_count * log(rule_count)

    # Calculate rule probabilities for potential use in diversity assessment
    probs <- sapply(grammar$rules[[nt]], function(rule) rule$prob)

    # For potential future extensions, store probability information
    attr(analysis$rule_diversity_index[[nt]], "probabilities") <- probs
  }

  # 5. Calculate Semantic-Aware RDI if requested
  if (use_semantic_awareness) {
    for (nt in grammar$non_terminals) {
      # Skip if non-terminal has no rules defined
      if (!(nt %in% names(grammar$rules)) || length(grammar$rules[[nt]]) <= 1) {
        analysis$semantic_rule_diversity_index[[nt]] <- 0
        next
      }

      # Group rules by semantic equivalence
      semantic_groups <- semantic_group_fn(grammar$rules[[nt]], nt)

      # Count the number of semantically distinct groups
      semantic_rule_count <- length(semantic_groups)

      if (semantic_rule_count <= 1) {
        analysis$semantic_rule_diversity_index[[nt]] <- 0
        next
      }

      # Calculate SARDI as semantic_rule_count * log(semantic_rule_count)
      analysis$semantic_rule_diversity_index[[nt]] <- semantic_rule_count * log(semantic_rule_count)
    }
  }

  # 6. Sort non-terminals by diversity potential
  # Create a data frame for easier sorting
  diversity_df <- data.frame(
    non_terminal = character(0),
    rule_count = numeric(0),
    diversity_index = numeric(0),
    stringsAsFactors = FALSE
  )

  # Populate data frame
  for (nt in names(analysis$rule_diversity_index)) {
    diversity_index <- if (use_semantic_awareness && !is.null(analysis$semantic_rule_diversity_index[[nt]])) {
      analysis$semantic_rule_diversity_index[[nt]]
    } else {
      analysis$rule_diversity_index[[nt]]
    }

    diversity_df <- rbind(diversity_df, data.frame(
      non_terminal = nt,
      rule_count = analysis$rule_counts[[nt]],
      diversity_index = diversity_index,
      stringsAsFactors = FALSE
    ))
  }

  # Sort by diversity index (descending)
  diversity_df <- diversity_df[order(-diversity_df$diversity_index), ]
  analysis$sorted_by_diversity <- diversity_df$non_terminal

  # Sort by rule count (descending)
  rule_count_df <- diversity_df[order(-diversity_df$rule_count), ]
  analysis$sorted_by_rules <- rule_count_df$non_terminal

  # 7. Identify high-diversity grammar paths
  # Start with the most diverse non-terminals
  top_diverse_nts <- head(analysis$sorted_by_diversity, min(5, length(analysis$sorted_by_diversity)))

  for (nt in top_diverse_nts) {
    if (!(nt %in% names(grammar$rules))) next

    # Create a list to store paths for this non-terminal
    analysis$high_diversity_paths[[nt]] <- list()

    # If this is the start symbol, include all rule indices
    if (nt == grammar$start_symbol) {
      analysis$high_diversity_paths[[nt]]$all_rules <- seq_along(grammar$rules[[nt]])
    } else {
      # For other high-diversity non-terminals, categorize rules

      # Get non-recursive rules using the package's existing function
      non_recursive_rules <- get_non_recursive_rules(grammar_with_flags, nt)
      if (length(non_recursive_rules) > 0) {
        # Find indices of non-recursive rules in the original grammar
        nr_indices <- integer(0)
        for (nr_rule in non_recursive_rules) {
          for (i in seq_along(grammar$rules[[nt]])) {
            if (identical(nr_rule$symbols, grammar$rules[[nt]][[i]]$symbols) &&
              identical(nr_rule$types, grammar$rules[[nt]][[i]]$types)) {
              nr_indices <- c(nr_indices, i)
              break
            }
          }
        }
        analysis$high_diversity_paths[[nt]]$non_recursive_rules <- nr_indices
      }

      # Find recursive rules (complement of non-recursive rules)
      all_indices <- seq_along(grammar$rules[[nt]])
      recursive_indices <- setdiff(all_indices, analysis$high_diversity_paths[[nt]]$non_recursive_rules)
      if (length(recursive_indices) > 0) {
        analysis$high_diversity_paths[[nt]]$recursive_rules <- recursive_indices
      }

      # Calculate rule complexity (simple heuristic: count of NT symbols)
      rule_complexity <- sapply(grammar$rules[[nt]], function(rule) {
        sum(rule$types == "NT")
      })

      # Select rules with above-average complexity
      above_avg <- which(rule_complexity > mean(rule_complexity))
      if (length(above_avg) > 0) {
        analysis$high_diversity_paths[[nt]]$complex_rules <- above_avg
      }

      # Also include diverse but simple rules
      below_avg <- which(rule_complexity <= mean(rule_complexity))
      if (length(below_avg) > 0) {
        analysis$high_diversity_paths[[nt]]$simple_rules <- below_avg
      }
    }
  }

  return(analysis)
}

#' Create Directed Genotype Templates from Grammar Analysis
#'
#' Produces a set of partial genotype "templates" that bias
#' the GE mapping process toward high-diversity grammar paths,
#' with particular emphasis on less common rules.
#'
#' @param analysis The list returned by analyze_grammar_for_diversity().
#' @param max_diverse_nt Integer: how many non-terminals to target (in
#'   order of descending diversity). Defaults to 5.
#' @param max_rules_per_nt Integer: maximum number of distinct rules
#'   to create templates for in each targeted non-terminal. Defaults to 3.
#' @param include_start_symbol Logical: whether to always create templates
#'   that select each rule of the start symbol. Defaults to TRUE.
#' @param prioritize_uncommon Logical: whether to prioritize less common rules
#'   (lower probability rules) when selecting templates. Defaults to TRUE.
#' @param create_combined_templates Logical: whether to create templates that
#'   combine choices for multiple non-terminals. Defaults to TRUE.
#' @param dedup_templates Logical: whether to remove duplicate templates.
#'   Defaults to TRUE.
#'
#' @return A list of partial templates. Each element is itself a named list,
#'   e.g. \code{list(expr = 0.3, product_expr = 0.6, ...)} indicating
#'   codons for some non-terminals. Others remain \code{NULL} (or omitted),
#'   allowing for random filling later.
#'
#' @export
directed_templates_from_analysis <- function(
    analysis,
    max_diverse_nt = 5,
    max_rules_per_nt = 3,
    include_start_symbol = TRUE,
    prioritize_uncommon = TRUE,
    create_combined_templates = TRUE,
    dedup_templates = TRUE
) {
  # Helper to compute a codon that forces selection of a specific rule
  pick_mean_codon <- function(prob_vec, rule_index) {
    # prob_vec: e.g., c(0.2, 0.2, 0.2, 0.2, 0.2)
    # rule_index: 1-based index of the target rule
    cum_probs <- cumsum(prob_vec)
    lower_bound <- if (rule_index == 1) 0 else cum_probs[rule_index - 1]
    upper_bound <- cum_probs[rule_index]
    # Return the midpoint of the probability range for this rule
    (lower_bound + upper_bound) / 2
  }

  # Helper to check if a template is a duplicate of existing ones
  is_duplicate <- function(new_tmpl, existing_templates) {
    for (tmpl in existing_templates) {
      # Different number of non-terminals means not a duplicate
      if (length(tmpl) != length(new_tmpl)) next
      # Check if all non-terminals and their values match
      if (all(names(tmpl) == names(new_tmpl)) &&
          all(sapply(names(tmpl), function(n) tmpl[[n]] == new_tmpl[[n]]))) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  # Grab top N non-terminals by diversity
  sorted_nts <- analysis$sorted_by_diversity
  targeted_nts <- head(sorted_nts, max_diverse_nt)

  # Always include the start symbol if requested (and not already in top N)
  if (isTRUE(include_start_symbol)) {
    start_nt <- analysis$start_symbol
    if (!start_nt %in% targeted_nts) {
      targeted_nts <- unique(c(start_nt, targeted_nts))
    }
  }

  all_templates <- list()

  # For each targeted non-terminal:
  for (nt in targeted_nts) {
    # Retrieve the associated probabilities from the analysis
    rdi_entry <- analysis$rule_diversity_index[[nt]]
    if (is.null(rdi_entry) || rdi_entry == 0) {
      # Means it has <= 1 rule or is otherwise uninteresting for templates
      next
    }
    prob_vec <- attr(rdi_entry, "probabilities")
    if (is.null(prob_vec) || length(prob_vec) < 2) {
      next
    }

    # Calculate relative diversity for this non-terminal
    # (Used to scale the number of templates we create)
    relative_diversity <- 1.0
    if (length(analysis$rule_diversity_index) > 0) {
      all_rdi_values <- unlist(analysis$rule_diversity_index)
      if (length(all_rdi_values) > 0 && max(all_rdi_values) > 0) {
        relative_diversity <- rdi_entry / max(all_rdi_values)
      }
    }
    # Scale the number of templates based on relative diversity (min 1)
    nt_max_rules <- max(1, round(max_rules_per_nt * relative_diversity))

    # Decide which specific rules to focus on
    rule_count <- length(prob_vec)
    candidate_rule_indices <- seq_len(rule_count)

    # If we have high_diversity_paths info for nt, use that to filter candidate rules
    chosen_rules <- candidate_rule_indices
    if (!is.null(analysis$high_diversity_paths[[nt]])) {
      hdp <- analysis$high_diversity_paths[[nt]]
      # Collect all known rule subsets
      subsets <- c(
        hdp$all_rules,
        hdp$non_recursive_rules,
        hdp$recursive_rules,
        hdp$complex_rules,
        hdp$simple_rules
      )
      subsets <- unique(subsets[!is.na(subsets)])
      # Only keep valid rule indices
      subsets <- subsets[subsets >= 1 & subsets <= rule_count]
      if (length(subsets) > 0) {
        chosen_rules <- subsets
      }
    }

    # Sophisticated rule selection based on prioritize_uncommon parameter
    if (length(chosen_rules) > nt_max_rules) {
      if (isTRUE(prioritize_uncommon)) {
        # Sort by probability (ascending) to prioritize less common rules
        # This increases diversity by ensuring rare paths are explored
        sorted_indices <- order(prob_vec[chosen_rules])
        chosen_rules <- chosen_rules[sorted_indices[seq_len(nt_max_rules)]]
      } else {
        # Evenly distribute choices across the range of rules
        # This ensures good coverage across the rule space
        chosen_rules <- chosen_rules[round(seq(1, length(chosen_rules), length.out = nt_max_rules))]
      }
    }

    # For each chosen rule, create a template that sets nt's codon to force that rule
    for (rule_i in chosen_rules) {
      codon <- pick_mean_codon(prob_vec, rule_i)
      # Build a partial template
      tmpl <- list()
      tmpl[[nt]] <- codon

      # Add template if it's not a duplicate (if deduplication is enabled)
      if (!isTRUE(dedup_templates) || !is_duplicate(tmpl, all_templates)) {
        all_templates <- c(all_templates, list(tmpl))
      }
    }
  }

  # Additionally, if we are including the start symbol and it has 'all_rules' in
  # high_diversity_paths, ensure each rule is represented at least once.
  if (include_start_symbol && !is.null(analysis$high_diversity_paths[[analysis$start_symbol]]$all_rules)) {
    start_symbol <- analysis$start_symbol
    path_info <- analysis$high_diversity_paths[[start_symbol]]
    rule_indices <- path_info$all_rules
    if (!is.null(rule_indices) && length(rule_indices) > 1) {
      start_prob_vec <- attr(analysis$rule_diversity_index[[start_symbol]], "probabilities")
      for (ri in rule_indices) {
        codon <- pick_mean_codon(start_prob_vec, ri)
        tmpl <- list()
        tmpl[[start_symbol]] <- codon

        # Add template if it's not a duplicate
        if (!isTRUE(dedup_templates) || !is_duplicate(tmpl, all_templates)) {
          all_templates <- c(all_templates, list(tmpl))
        }
      }
    }
  }

  # Create combined templates for pairs of high-diversity non-terminals
  if (isTRUE(create_combined_templates) && length(targeted_nts) >= 2) {
    # Limit the number of combined templates to avoid explosion
    max_combined_templates <- min(10, length(targeted_nts) * (length(targeted_nts) - 1) / 2)
    combined_count <- 0

    # Try pairs of non-terminals in order of diversity
    for (i in 1:(length(targeted_nts) - 1)) {
      if (combined_count >= max_combined_templates) break

      nt1 <- targeted_nts[i]
      for (j in (i+1):length(targeted_nts)) {
        if (combined_count >= max_combined_templates) break

        nt2 <- targeted_nts[j]

        # Get probability vectors for both non-terminals
        prob_vec1 <- attr(analysis$rule_diversity_index[[nt1]], "probabilities")
        prob_vec2 <- attr(analysis$rule_diversity_index[[nt2]], "probabilities")

        if (!is.null(prob_vec1) && !is.null(prob_vec2) &&
            length(prob_vec1) > 1 && length(prob_vec2) > 1) {

          # For combined templates, prioritize unusual combinations
          # by choosing low-probability rules
          if (isTRUE(prioritize_uncommon)) {
            # Select rules weighted inversely to their probability
            # (makes rare rules more likely to be selected)
            weights1 <- 1 / (prob_vec1 + 0.01)  # Add small constant to avoid division by zero
            weights2 <- 1 / (prob_vec2 + 0.01)

            rule1 <- sample(seq_along(prob_vec1), 1, prob = weights1 / sum(weights1))
            rule2 <- sample(seq_along(prob_vec2), 1, prob = weights2 / sum(weights2))
          } else {
            # Uniform random selection
            rule1 <- sample(seq_along(prob_vec1), 1)
            rule2 <- sample(seq_along(prob_vec2), 1)
          }

          # Create the combined template
          tmpl <- list()
          tmpl[[nt1]] <- pick_mean_codon(prob_vec1, rule1)
          tmpl[[nt2]] <- pick_mean_codon(prob_vec2, rule2)

          # Add template if it's not a duplicate
          if (!isTRUE(dedup_templates) || !is_duplicate(tmpl, all_templates)) {
            all_templates <- c(all_templates, list(tmpl))
            combined_count <- combined_count + 1
          }
        }
      }
    }
  }

  # Return list of partial templates
  return(all_templates)
}

#' Create a Single Individual from a Directed Template
#'
#' Builds a genotype by combining the user-supplied codons in \code{template}
#' with random codons for any other non-terminals, then maps the genotype to
#' a phenotype. Returns a standard "individual" list.
#'
#' @param grammar A PCFG grammar structure (e.g., from \code{read_pcfg_grammar()}).
#' @param template A named list of codons (e.g. \code{list(expr = 0.25, product_expr = 0.6)})
#'   where each codon is numeric and indicates how to select a particular rule.
#' @param max_depth Integer maximum derivation depth. Passed to \code{map_genotype()}.
#' @param random_codons_per_nt Integer: how many codons to randomly generate
#'   for non-terminals not in the template. Default 1.
#' @param lazy_codons Logical: whether to create random codons only when needed during
#'   mapping (TRUE) or pre-generate them all (FALSE). Defaults to FALSE.
#' @param seed Integer: optional random seed for reproducibility.
#' @return A list representing the new individual, including genotype, phenotype, etc.
#'
#' @export
create_individual_from_template <- function(grammar,
                                            template,
                                            max_depth = 17,
                                            random_codons_per_nt = 1,
                                            lazy_codons = FALSE,
                                            seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit(.Random.seed <- old_seed)
    set.seed(seed)
  }

  # For lazy codon generation, we only initialize the genotype with template values
  # and let map_genotype generate other values as needed
  if (isTRUE(lazy_codons)) {
    genotype <- list()
    for (nt in names(template)) {
      codon_value <- template[[nt]]
      genotype[[nt]] <- if (length(codon_value) == 1) {
        c(codon_value)
      } else {
        as.numeric(codon_value)
      }
    }
  } else {
    # Prepare a complete genotype with values for all non-terminals
    nt_list <- grammar$non_terminals
    genotype <- vector("list", length(nt_list))
    names(genotype) <- nt_list

    # Fill in codons from template if present, otherwise generate random codons
    for (nt in nt_list) {
      if (nt %in% names(template)) {
        # Use codon(s) from template
        codon_value <- template[[nt]]
        genotype[[nt]] <- if (length(codon_value) == 1) {
          c(codon_value)
        } else {
          as.numeric(codon_value)
        }
      } else {
        # Generate random codons in [0,1]
        genotype[[nt]] <- runif(random_codons_per_nt, 0, 1)
      }
    }
  }

  # Map genotype to phenotype
  # Note: If lazy_codons is TRUE, map_genotype must be able to handle generating
  # random codons for non-terminals not in the genotype
  mapping_result <- map_genotype(grammar, genotype, max_depth)

  # Build the final "individual" object
  individual <- list(
    genotype = genotype,
    phenotype = mapping_result$phenotype,
    fitness = NULL,  # Not yet evaluated
    mapping_positions = mapping_result$positions,
    tree_depth = mapping_result$depth
  )

  return(individual)
}

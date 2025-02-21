#' Update grammar probabilities based on a successful individual
#'
#' @param grammar PCFG grammar structure
#' @param individual Individual with genotype and positions
#' @param learning_factor Learning rate for probability adjustment (default: 0.01)
#' @param force_perturbation Whether to force perturbation of all probability distributions (default: FALSE)
#' @param perturbation_factor Amount of perturbation when forcing changes (default: 0.1)
#' @return Updated grammar with adjusted probabilities
#' @export
update_grammar_probabilities <- function(grammar, individual, learning_factor = 0.01,
                                         force_perturbation = FALSE, perturbation_factor = 0.1) {
  if (is.null(individual$genotype)) {
    stop("Invalid individual: missing genotype")
  }

  if (length(individual$genotype) == 0) {
    stop("Invalid individual: empty genotype")
  }

  # Create a deep copy with explicit new lists to avoid reference issues
  grammar_copy <- list(
    non_terminals = grammar$non_terminals,
    terminals = grammar$terminals,
    start_symbol = grammar$start_symbol,
    rules = list()
  )

  # Deep copy of rules
  for (nt in names(grammar$rules)) {
    grammar_copy$rules[[nt]] <- list()
    for (i in seq_along(grammar$rules[[nt]])) {
      grammar_copy$rules[[nt]][[i]] <- list(
        symbols = grammar$rules[[nt]][[i]]$symbols,
        types = grammar$rules[[nt]][[i]]$types,
        prob = grammar$rules[[nt]][[i]]$prob
      )
      if (!is.null(grammar$rules[[nt]][[i]]$is_recursive)) {
        grammar_copy$rules[[nt]][[i]]$is_recursive <- grammar$rules[[nt]][[i]]$is_recursive
      }
    }
  }

  # Process every non-terminal with a fixed perturbation if requested
  if (force_perturbation) {
    for (nt in names(grammar_copy$rules)) {
      if (length(grammar_copy$rules[[nt]]) > 1) {
        # Find the highest probability rule
        highest_prob <- 0
        highest_idx <- 1

        for (i in seq_along(grammar_copy$rules[[nt]])) {
          if (grammar_copy$rules[[nt]][[i]]$prob > highest_prob) {
            highest_prob <- grammar_copy$rules[[nt]][[i]]$prob
            highest_idx <- i
          }
        }

        # Ensure a meaningful change
        perturb_amount <- min(perturbation_factor, highest_prob * 0.5)

        # Reduce the highest probability
        grammar_copy$rules[[nt]][[highest_idx]]$prob <- grammar_copy$rules[[nt]][[highest_idx]]$prob - perturb_amount

        # Find a random rule to increase (not the highest)
        other_indices <- setdiff(seq_along(grammar_copy$rules[[nt]]), highest_idx)
        if (length(other_indices) > 0) {
          lucky_idx <- sample(other_indices, 1)
          grammar_copy$rules[[nt]][[lucky_idx]]$prob <- grammar_copy$rules[[nt]][[lucky_idx]]$prob + perturb_amount
        } else {
          # If there's only one rule, restore its probability
          grammar_copy$rules[[nt]][[highest_idx]]$prob <- 1.0
        }

        # Normalize again
        total_prob <- sum(sapply(grammar_copy$rules[[nt]], function(r) r$prob))
        for (i in seq_along(grammar_copy$rules[[nt]])) {
          grammar_copy$rules[[nt]][[i]]$prob <- grammar_copy$rules[[nt]][[i]]$prob / total_prob
        }
      }
    }
  }

  # Now perform the regular probability updates for used non-terminals
  used_nts <- names(individual$genotype)[sapply(individual$genotype, length) > 0]
  for (nt in used_nts) {
    if (!(nt %in% names(grammar_copy$rules))) {
      next
    }

    # Count rule usage
    rule_counts <- integer(length(grammar_copy$rules[[nt]]))
    for (codon in individual$genotype[[nt]]) {
      cum_prob <- 0
      for (i in seq_along(grammar_copy$rules[[nt]])) {
        cum_prob <- cum_prob + grammar_copy$rules[[nt]][[i]]$prob
        if (codon <= cum_prob) {
          rule_counts[i] <- rule_counts[i] + 1
          break
        }
      }
    }

    # Update probabilities further based on usage
    total_used <- sum(rule_counts)
    if (total_used > 0) {
      # Apply learning factor adjustments
      for (i in seq_along(grammar_copy$rules[[nt]])) {
        if (rule_counts[i] > 0) {
          grammar_copy$rules[[nt]][[i]]$prob <- grammar_copy$rules[[nt]][[i]]$prob +
            (learning_factor * (rule_counts[i] / total_used))
        } else {
          grammar_copy$rules[[nt]][[i]]$prob <- grammar_copy$rules[[nt]][[i]]$prob *
            (1 - learning_factor)
        }
      }

      # Normalize again
      total_prob <- sum(sapply(grammar_copy$rules[[nt]], function(r) r$prob))
      for (i in seq_along(grammar_copy$rules[[nt]])) {
        grammar_copy$rules[[nt]][[i]]$prob <- grammar_copy$rules[[nt]][[i]]$prob / total_prob
      }
    }
  }

  return(grammar_copy)
}

#' Initialize a population of random individuals
#'
#' @param grammar PCFG grammar structure
#' @param pop_size Population size
#' @param max_depth Maximum tree depth
#' @return List of individuals with genotypes
#' @export
initialize_population <- function(grammar, pop_size, max_depth = 17) {
  population <- vector("list", pop_size)

  for (i in seq_len(pop_size)) {
    population[[i]] <- list(
      genotype = generate_random_individual(grammar, max_depth),
      fitness = NULL
    )
  }

  return(population)
}

#' Evaluate all individuals in a population
#'
#' @param population List of individuals
#' @param grammar PCFG grammar structure
#' @param eval_func Fitness evaluation function
#' @param max_depth Maximum tree depth
#' @param parallel Boolean indicating whether to use parallel evaluation
#' @param n_workers Number of workers for parallel processing (default: 2)
#' @param future_strategy Future strategy to use (default: "multisession")
#'
#' @details
#' Note that we use progressr progress bars with this function:
#'
#' ```r
#' library(progressr)
#' handlers(global = TRUE)
#' handlers("progress")
#'
#' with_progress({
#'   population <- evaluate_population(
#'     population = initial_pop,
#'     grammar = my_grammar,
#'     eval_func = my_fitness_function,
#'     parallel = TRUE
#'   )
#' })
#' ```
#'
#' @return Evaluated population with fitness values
#' @export
evaluate_population <- function(population, grammar, eval_func, max_depth = 17,
                                parallel = FALSE, n_workers = 2,
                                future_strategy = "multisession") {

  # Create progress handler for this function
  p <- progressr::progressor(steps = length(population))

  if (parallel) {
    # Setup parallel backend
    future::plan(future_strategy, workers = n_workers)
    doFuture::registerDoFuture()

    # Use foreach for parallel evaluation
    evaluated_population <- foreach::foreach(
      ind = population,
      .packages = c("Rpsge"), # Add any other required packages
      .export = c("map_genotype", "eval_func") # Add other functions that need to be exported
      ) %dorng% {
      ind_copy <- ind

      # Map genotype to phenotype
      mapping_result <- map_genotype(grammar, ind_copy$genotype, max_depth)
      ind_copy$phenotype <- mapping_result$phenotype
      ind_copy$mapping_positions <- mapping_result$positions
      ind_copy$tree_depth <- mapping_result$depth

      # Evaluate fitness
      ind_copy$fitness <- eval_func(ind_copy$phenotype)

      # Update progress (works in parallel with progressr)
      p(sprintf("Evaluated individual: fitness = %.6f", ind_copy$fitness))

      return(ind_copy)
    }

    # Convert to list if not already
    evaluated_population <- as.list(evaluated_population)

  } else {
    # Sequential evaluation
    evaluated_population <- lapply(population, function(ind) {
      ind_copy <- ind

      # Map genotype to phenotype
      mapping_result <- map_genotype(grammar, ind_copy$genotype, max_depth)
      ind_copy$phenotype <- mapping_result$phenotype
      ind_copy$mapping_positions <- mapping_result$positions
      ind_copy$tree_depth <- mapping_result$depth

      # Evaluate fitness
      ind_copy$fitness <- eval_func(ind_copy$phenotype)

      # Update progress
      p(sprintf("Evaluated individual: fitness = %.6f", ind_copy$fitness))

      return(ind_copy)
    })
  }

  return(evaluated_population)
}

#' Tournament selection for PSGE
#'
#' @param population List of evaluated individuals
#' @param tournament_size Number of individuals to include in each tournament
#' @return Selected individual
#' @export
tournament_selection <- function(population, tournament_size = 3) {
  if (tournament_size > length(population)) {
    tournament_size <- length(population)
  }

  # Randomly select tournament_size individuals
  candidates_indices <- sample(seq_along(population), tournament_size, replace = FALSE)
  candidates <- population[candidates_indices]

  # Find the individual with best fitness (minimum value)
  fitness_values <- sapply(candidates, function(ind) ind$fitness)
  best_index <- which.min(fitness_values)

  # Return a copy of the selected individual to prevent modification of the original
  selected <- candidates[[best_index]]
  return(selected)
}

#' Crossover operator for PSGE
#'
#' @param parent1 First parent individual
#' @param parent2 Second parent individual
#' @return New individual created by crossover
#' @export
crossover <- function(parent1, parent2) {
  # Create binary mask for each non-terminal
  all_nts <- unique(c(names(parent1$genotype), names(parent2$genotype)))
  mask <- setNames(
    sample(0:1, length(all_nts), replace = TRUE),
    all_nts
  )

  # Initialize new genotype
  offspring_genotype <- list()

  # Apply crossover based on mask
  for (nt in all_nts) {
    if (mask[nt] == 1) {
      if (nt %in% names(parent1$genotype)) {
        offspring_genotype[[nt]] <- parent1$genotype[[nt]]
      } else {
        offspring_genotype[[nt]] <- numeric(0)
      }
    } else {
      if (nt %in% names(parent2$genotype)) {
        offspring_genotype[[nt]] <- parent2$genotype[[nt]]
      } else {
        offspring_genotype[[nt]] <- numeric(0)
      }
    }
  }

  # Create offspring with NULL fitness to ensure re-evaluation
  offspring <- list(
    genotype = offspring_genotype,
    fitness = NULL
  )

  return(offspring)
}

#' Mutation operator for PSGE
#'
#' @param individual Individual to mutate
#' @param mutation_rate Probability of mutating each codon
#' @return Mutated individual
#' @export
mutate <- function(individual, mutation_rate = 0.1) {
  # Create a deep copy of the genotype
  new_genotype <- lapply(individual$genotype, function(values) {
    if (length(values) > 0) {
      return(values)
    } else {
      return(numeric(0))
    }
  })

  # Flag to track if any mutation occurred
  any_mutation <- FALSE

  # Apply mutation to each non-terminal
  for (nt in names(new_genotype)) {
    for (i in seq_along(new_genotype[[nt]])) {
      if (runif(1) < mutation_rate) {
        # Apply Gaussian mutation
        mutation <- rnorm(1, mean = 0, sd = 0.5)
        new_value <- new_genotype[[nt]][i] + mutation

        # Keep value in [0,1] range
        new_genotype[[nt]][i] <- max(0, min(1, new_value))

        any_mutation <- TRUE
      }
    }
  }

  # Create new individual with NULL fitness to ensure re-evaluation
  mutated <- list(
    genotype = new_genotype,
    fitness = NULL
  )

  return(mutated)
}

#' Create the next generation through selection, crossover, and mutation
#'
#' @param population Current population
#' @param grammar PCFG grammar structure
#' @param crossover_rate Probability of applying crossover
#' @param mutation_rate Probability of mutating each codon
#' @param tournament_size Tournament selection size
#' @param elitism Number of top individuals to preserve
#' @return New population
#' @export
create_next_generation <- function(population, grammar, crossover_rate = 0.9,
                                   mutation_rate = 0.1, tournament_size = 3, elitism = 10) {
  # Sort population by fitness
  sorted_pop <- population[order(sapply(population, function(ind) ind$fitness))]

  # Determine elite count (handle edge case where elitism exceeds population size)
  elite_count <- min(elitism, length(sorted_pop))

  # Create output population vector
  new_pop_size <- length(sorted_pop)
  new_population <- vector("list", new_pop_size)

  # First create the non-elite individuals
  for (i in 1:(new_pop_size - elite_count)) {
    if (runif(1) < crossover_rate) {
      # Crossover
      parent1 <- tournament_selection(sorted_pop, tournament_size)
      parent2 <- tournament_selection(sorted_pop, tournament_size)
      offspring <- crossover(parent1, parent2)
    } else {
      # Copy but reset fitness
      selected <- tournament_selection(sorted_pop, tournament_size)
      offspring <- list(
        genotype = selected$genotype,
        fitness = NULL
      )
    }

    # Apply mutation
    offspring <- mutate(offspring, mutation_rate)

    # Ensure fitness is NULL
    offspring$fitness <- NULL

    # Add to new population
    new_population[[i]] <- offspring
  }

  # Now add the elite individuals at the end
  if (elite_count > 0) {
    for (i in 1:elite_count) {
      elite_index <- i  # Get from the beginning (best individuals)
      new_index <- new_pop_size - elite_count + i  # Place at the end

      # Create a copy of the elite individual
      elite_copy <- sorted_pop[[elite_index]]

      # Place in new population
      new_population[[new_index]] <- elite_copy
    }
  }

  return(new_population)
}

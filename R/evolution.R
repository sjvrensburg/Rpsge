# Helper functions
sort_by_fitness <- function(pop) {
  pop[order(sapply(pop, `[[`, "fitness"))]
}

Evolution <- R6::R6Class(
  "Evolution",
  public = list(
    grammar = NULL,
    eval_func = NULL,
    params = list(
      popsize = 100,
      generations = 50,
      elitism = 10,
      prob_crossover = 0.9,
      prob_mutation = 0.1,
      tournament_size = 3,
      learning_factor = 0.01,
      adaptive_lf = FALSE,
      adaptive_increment = 0.0001,
      seed = NULL
    ),

    initialize = function(grammar, eval_func, params = list()) {
      self$grammar <- grammar
      self$eval_func <- eval_func

      # Override default parameters
      for (param_name in names(params)) {
        self$params[[param_name]] <- params[[param_name]]
      }

      if (!is.null(self$params$seed)) {
        set.seed(self$params$seed)
      }
    },

    generate_random_individual = function() {
      genotype <- vector("list", length = length(self$grammar$non_terminals$values()))
      names(genotype) <- self$grammar$non_terminals$values()

      result <- self$grammar$recursive_individual_creation(
        genotype = genotype,
        symbol = self$grammar$start_rule[[1]],
        current_depth = 0
      )

      list(
        genotype = result$genotype,
        fitness = NULL,
        tree_depth = result$depth
      )
    },

    evaluate = function(ind) {
      # Initialize mapping positions
      positions_to_map <- rep(0, length(ind$genotype))

      # Map genotype to phenotype
      mapping_result <- self$grammar$mapping(ind$genotype, positions_to_map)

      # Evaluate fitness using provided function
      quality <- self$eval_func(mapping_result[[1]])  # First element is phenotype

      # Update individual with results
      ind$phenotype <- mapping_result[[1]]
      ind$fitness <- quality
      ind$mapping_values <- positions_to_map

      ind
    },

    tournament_selection = function(population) {
      # Randomly select tournament_size individuals
      candidates <- sample(population, self$params$tournament_size, replace = FALSE)

      # Return the individual with best fitness (minimum value)
      candidates[order(sapply(candidates, function(x) x$fitness))][[1]]
    },

    validate_phenotype = function(genotype) {
      positions_to_map <- rep(0, length(genotype))
      tryCatch({
        self$grammar$mapping(genotype, positions_to_map)
        TRUE
      }, error = function(e) FALSE)
    },

    crossover_impl = function(p1, p2) {
      # Initialize new genotype
      genotype <- list()

      # For each non-terminal in grammar
      for (nt in names(p1$genotype)) {
        current_nt <- list()
        nt_index <- self$grammar$index_of_non_terminal[[nt]]

        # Get production counts for both parents
        p1_productions <- sapply(p1$genotype[[nt]], function(x) x[1])
        p2_productions <- sapply(p2$genotype[[nt]], function(x) x[1])

        # Get PCFG probabilities
        probs <- self$grammar$pcfg[nt_index,]

        # For each position that exists in either parent
        max_len <- max(length(p1_productions), length(p2_productions))
        for (i in seq_len(max_len)) {
          if (i <= length(p1_productions) && i <= length(p2_productions)) {
            # Both parents have this position
            p1_prob <- probs[p1_productions[i]]
            p2_prob <- probs[p2_productions[i]]

            # Select based on relative probabilities
            total_prob <- p1_prob + p2_prob
            if (runif(1) < p1_prob/total_prob) {
              current_nt[[i]] <- p1$genotype[[nt]][[i]]
            } else {
              current_nt[[i]] <- p2$genotype[[nt]][[i]]
            }
          } else if (i <= length(p1_productions)) {
            # Only p1 has this position
            current_nt[[i]] <- p1$genotype[[nt]][[i]]
          } else {
            # Only p2 has this position
            current_nt[[i]] <- p2$genotype[[nt]][[i]]
          }
        }

        genotype[[nt]] <- current_nt
      }

      names(genotype) <- names(p1$genotype)

      list(
        genotype = genotype,
        fitness = NA,
        mapping_values = rep(0, length(genotype))
      )
    },

    crossover = function(p1, p2) {
      attempts <- 0
      repeat {
        offspring <- self$crossover_impl(p1, p2)  # Current crossover logic
        if (self$validate_phenotype(offspring$genotype) || attempts > 10) break
        attempts <- attempts + 1
      }
      offspring
    },

    mutate_impl = function(ind) {
      # Mutate each non-terminal's rules
      for (nt in names(ind$genotype)) {
        for (i in seq_along(ind$genotype[[nt]])) {
          if (runif(1) < self$params$prob_mutation) {
            # Get current values
            current_value <- ind$genotype[[nt]][[i]]
            nt_index <- self$grammar$index_of_non_terminal[[nt]]

            # Get probability distribution for this non-terminal
            probs <- self$grammar$pcfg[nt_index,]

            # Calculate variance based on probability distribution
            # Using entropy as a measure of uncertainty
            entropy <- -sum(probs * log(probs + 1e-10))
            variance <- min(0.5, entropy / 2)  # Cap at 0.5

            # Apply Gaussian mutation with dynamic variance
            codon <- stats::rnorm(1, current_value[2], variance)
            codon <- max(0, min(1, codon))

            # Select rule based on mutated codon
            prob_aux <- 0
            for (index in seq_along(self$grammar$grammar[[nt]])) {
              prob_aux <- prob_aux + self$grammar$pcfg[nt_index, index]
              if (codon <= round(prob_aux, 3)) {
                ind$genotype[[nt]][[i]] <- c(index, codon, current_value[3])
                break
              }
            }
          }
        }
      }

      ind$fitness <- NA
      ind
    },

    mutate = function(p) {
      attempts <- 0
      repeat {
        offspring <- self$mutate_impl(p)  # Current mutate logic
        if (self$validate_phenotype(offspring$genotype) || attempts > 10) break
        attempts <- attempts + 1
      }
      offspring
    },

    update_pcfg = function(best_individual, learning_factor) {
      # For each non-terminal in the grammar
      for (nt in names(best_individual$genotype)) {
        # Get the index for this non-terminal
        nt_index <- self$grammar$index_of_non_terminal[[nt]]

        # Get the number of possible production rules for this non-terminal
        n_rules <- ncol(self$grammar$pcfg)

        # Initialize counter vector for rules
        counter <- integer(n_rules)

        # Count occurrences of each production rule
        for (gene in best_individual$genotype[[nt]]) {
          rule_index <- gene[1]
          counter[rule_index] <- counter[rule_index] + 1
        }

        # Calculate total number of rules used
        total <- sum(counter)

        if (total > 0) {
          # Get current probabilities for this non-terminal
          current_probs <- self$grammar$pcfg[nt_index, ]

          # Update probabilities based on rule usage
          new_probs <- current_probs

          for (j in seq_along(counter)) {
            if (counter[j] > 0) {
              # Increase probability for used rules
              new_probs[j] <- min(current_probs[j] +
                                    learning_factor * (counter[j] / total), 1.0)
            } else {
              # Decrease probability for unused rules
              new_probs[j] <- max(current_probs[j] -
                                    learning_factor * current_probs[j], 0.0)
            }
          }

          # Normalize probabilities to ensure they sum to 1
          self$grammar$pcfg[nt_index, ] <- new_probs / sum(new_probs)
        }
      }
    },

    make_initial_population = function() {
      population <- vector("list", length = self$params$popsize)

      for (i in seq_len(self$params$popsize)) {
        individual <- self$generate_random_individual()
        population[[i]] <- individual
      }

      return(population)
    },

    run_evolution = function() {
      # Initialize and evaluate population
      message("Generating the initial population...")
      population <- self$make_initial_population()
      message("Evaluating the initial population...")
      population <- self$evaluate_all(population)
      message("Generated and evaluated the initial population...")
      # Tracking variables
      best_overall <- NULL
      best_generation <- NULL
      flag <- FALSE  # Flag for alternating PCFG updates

      # Main evolutionary loop
      message("Starting evolution...")
      for (gen in 1:self$params$generations) {
        message(sprintf("\tGen %d...", gen))
        message("\tSorting...")
        population <- sort_by_fitness(population)
        message("\tSorted...")
        # Update best individuals
        current_best <- population[[1]]
        message("\tCurrent best is:\n")
        print(current_best$phenotype)

        if (is.null(best_overall) || current_best$fitness < best_overall$fitness) {
          best_overall <- current_best
        }

        message("\tBest overal:\n")
        print(current_best$phenotype)

        # Update PCFG probabilities
        best_to_update <- if (!flag) best_overall else current_best
        message("\tUpdating PCFG probabilities")
        self$update_pcfg(best_to_update, self$params$learning_factor)
        message("\tUpdated PCFG probabilities")

        # Adaptive learning factor
        if (self$params$adaptive_lf) {
          self$params$learning_factor <-
            self$params$learning_factor + self$params$adaptive_increment
        }

        # Progress output
        cat(sprintf("Generation %d: Best Fitness = %f\n", gen, current_best$fitness))

        # Create new population
        new_population <- list()

        # Generate offspring
        while (length(new_population) < (self$params$popsize - self$params$elitism)) {
          # Select reproduction method
          if (runif(1) < self$params$prob_crossover) {
            # Crossover
            p1 <- self$tournament_selection(population)
            p2 <- self$tournament_selection(population)
            offspring <- self$crossover(p1, p2)
          } else {
            # Copy
            offspring <- self$tournament_selection(population)
          }

          # Apply mutation
          offspring <- self$mutate(offspring)

          # Evaluate offspring
          offspring <- self$evaluate(offspring)

          # Add to new population
          new_population <- c(new_population, list(offspring))
        }

        # Add elites to new population
        elites <- population[1:self$params$elitism]
        population <- c(new_population, elites)

        # Update tracking
        best_generation <- population[[1]]
        flag <- !flag  # Toggle flag for next generation
      }

      # Return best solution found
      return(best_overall)
    },
    evaluate_all = function(pop) {
      lapply(pop, self$evaluate)
    }
  )
)

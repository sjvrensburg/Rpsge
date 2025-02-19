# Helper functions
sort_by_fitness <- function(pop) {
  pop[order(sapply(pop, `[[`, "fitness"))]
}

Evolution <- R6::R6Class(
  "Evolution",
  private = list(
    validate_genotype = function(genotype) {
      if (!is.list(genotype)) {
        cat("Not a list\n")
        return(FALSE)
      }

      required_nts <- sapply(self$grammar$ordered_non_terminals$values(), `[[`, 1)
      if (!all(required_nts %in% names(genotype))) {
        cat("Missing NTs:", setdiff(required_nts, names(genotype)), "\n")
        return(FALSE)
      }

      for (nt in names(genotype)) {
        entries <- genotype[[nt]]
        if (!is.list(entries)) {
          cat("Entry for", nt, "not a list\n")
          return(FALSE)
        }
        for (entry in entries) {
          if (!is.numeric(entry) || length(entry) != 3) {
            cat("Invalid entry for", nt, ":", str(entry), "\n")
            return(FALSE)
          }
        }
      }
      TRUE
    }
  ),
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
    map_genotype = function(ind) {
      if (!private$validate_genotype(ind$genotype)) {
        return(NULL)
      }

      # Initialize positions
      positions_to_map <- rep(0, length(self$grammar$ordered_non_terminals$values()))

      # Map genotype using Grammar class methods
      result <- self$grammar$mapping(ind$genotype, positions_to_map)

      # Update individual with results
      ind$phenotype <- result$phenotype
      ind$mapping_values <- result$positions
      ind$tree_depth <- result$mapping_result

      ind
    },
    generate_random_individual = function() {
      # Initialize genotype with empty lists for all non-terminals
      genotype <- list()
      for (nt in names(self$grammar$grammar)) {
        genotype[[nt]] <- list()
      }

      # Create initial genotype starting from root
      result <- self$grammar$recursive_individual_creation(
        genotype = genotype,
        symbol = self$grammar$start_rule[[1]],
        current_depth = 0
      )

      # Return complete individual
      list(
        genotype = result$genotype,
        fitness = NULL,
        tree_depth = result$depth
      )
    },
    evaluate = function(ind) {
      ind <- self$map_genotype(ind)
      if (is.null(ind)) {
        stop("Invalid genotype during evaluation")
      }

      quality <- self$eval_func(ind$phenotype)
      ind$fitness <- quality
      ind$other_info <- list() # Maintain compatibility

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
      tryCatch(
        {
          self$grammar$mapping(genotype, positions_to_map)
          TRUE
        },
        error = function(e) FALSE
      )
    },
    crossover_impl = function(p1, p2) {
      # Create binary mask for each non-terminal
      mask <- sapply(names(p1$genotype), function(x) sample(0:1, 1))

      # Initialize new genotype
      genotype <- list()

      # For each non-terminal in grammar
      for (nt in names(p1$genotype)) {
        # Select parent based on mask
        genotype[[nt]] <- if (mask[nt] == 1) {
          p1$genotype[[nt]]
        } else {
          p2$genotype[[nt]]
        }
      }

      list(
        genotype = genotype,
        fitness = NA,
        mapping_values = rep(0, length(genotype))
      )
    },
    crossover = function(p1, p2) {
      attempts <- 0
      repeat {
        offspring <- self$crossover_impl(p1, p2) # Current crossover logic
        if (self$validate_phenotype(offspring$genotype) || attempts > 10) break
        attempts <- attempts + 1
      }
      offspring
    },
    mutate_impl = function(ind) {
      # For each non-terminal
      for (nt in names(ind$genotype)) {
        # For each codon in the non-terminal
        for (i in seq_along(ind$genotype[[nt]])) {
          if (runif(1) < self$params$prob_mutation) {
            # Get current values
            current_value <- ind$genotype[[nt]][[i]]

            # Apply Gaussian mutation with N(0, 0.50)
            codon <- current_value[2] + rnorm(1, 0, 0.50)
            codon <- max(0, min(1, codon)) # Keep in [0,1]

            # Calculate which rule this probability selects
            prob_aux <- 0
            nt_index <- self$grammar$index_of_non_terminal[[nt]]
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
        offspring <- self$mutate_impl(p) # Current mutate logic
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
      population <- self$evaluate_all(self$make_initial_population())
      # Tracking variables
      best_overall <- NULL
      best_generation <- NULL
      flag <- FALSE # Flag for alternating PCFG updates

      # Main evolutionary loop
      for (gen in 1:self$params$generations) {
        population <- sort_by_fitness(population)
        # Update best individuals
        current_best <- population[[1]]

        if (is.null(best_overall) || current_best$fitness < best_overall$fitness) {
          best_overall <- current_best
        }

        # Update PCFG probabilities
        best_to_update <- if (!flag) best_overall else current_best
        self$update_pcfg(best_to_update, self$params$learning_factor)

        # Adaptive learning factor
        if (self$params$adaptive_lf) {
          self$params$learning_factor <-
            self$params$learning_factor + self$params$adaptive_increment
        }

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
        flag <- !flag # Toggle flag for next generation
      }

      # Return best solution found
      return(best_overall)
    },
    evaluate_all = function(pop) {
      lapply(pop, self$evaluate)
    },
    print_genotype = function(ind) {
      str(ind$genotype)
    }
  )
)

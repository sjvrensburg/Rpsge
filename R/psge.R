#' Probabilistic Structured Grammatical Evolution
#'
#' @description
#' An R6 class implementing Probabilistic Structured Grammatical Evolution (PSGE).
#' PSGE combines Structured Grammatical Evolution with probabilistic grammar-guided
#' genetic programming to evolve solutions to various problems.
#'
#' @details
#' PSGE uses a probabilistic context-free grammar (PCFG) to guide the evolution
#' of solutions. The probabilities in the grammar are updated based on successful
#' individuals, biasing the search towards promising regions of the solution space.
#'
#' @examples
#' \dontrun{
#' # Create grammar
#' grammar <- create_grammar("path/to/grammar.bnf")
#'
#' # Define evaluation function
#' eval_func <- function(phenotype) {
#'   # Calculate and return fitness
#'   ...
#' }
#'
#' # Create PSGE instance
#' psge <- PSGE$new(
#'   grammar = grammar,
#'   eval_func = eval_func,
#'   population_size = 100
#' )
#'
#' # Run evolution
#' result <- psge$evolve()
#' }
#'
#' @export
PSGE <- R6::R6Class(
  "PSGE",

  private = list(
    .grammar = NULL,
    .population = NULL,
    .population_size = 100,
    .generations = 50,
    .elite_size = 10,
    .tournament_size = 3,
    .crossover_prob = 0.9,
    .mutation_prob = 0.1,
    .learning_factor = 0.01,
    .best_individual = NULL,
    .best_fitness = Inf,
    .evaluation_function = NULL,

    # Private methods are documented similar to public ones
    #' @description
    #' Generate a random individual
    #'
    #' @return List representing a new individual
    .generate_random_individual = function() {
      # Create empty genotype structure
      genotype <- replicate(
        length(private$.grammar$non_terminals),
        list(),
        simplify = FALSE
      )

      # Generate initial tree structure
      tree_depth <- private$.grammar$recursive_individual_creation(
        genome = genotype,
        symbol = private$.grammar$non_terminals[1],
        current_depth = 0
      )

      # Return structured individual
      list(
        genotype = genotype,
        fitness = NULL,
        tree_depth = tree_depth,
        mapping_values = numeric(length(private$.grammar$non_terminals))
      )
    },

    # Tournament selection
    .tournament_select = function() {
      # Handle case where tournament size > population size
      actual_tsize <- min(private$.tournament_size, length(private$.population))

      # Sample tournament participants
      pool <- sample(private$.population, actual_tsize)

      # Deep copy the winner
      winner <- pool[[which.min(sapply(pool, function(x) x$fitness))]]
      utils::modifyList(list(), winner)
    },

    # Gaussian mutation operator
    .mutate = function(individual, mutation_prob) {
      # Deep copy the individual
      ind <- utils::modifyList(list(), individual)
      ind$fitness <- NULL

      # Get grammar information
      pcfg <- private$.grammar$pcfg
      non_terminals <- private$.grammar$non_terminals

      # Find mutable genes (non-terminals with multiple production rules)
      for(nt_idx in seq_along(non_terminals)) {
        # Skip if no productions or only one rule
        if(length(ind$genotype[[nt_idx]]) == 0) next
        if(sum(pcfg[nt_idx,] > 0) <= 1) next

        # Check each position up to mapping value
        mapped <- ind$mapping_values[nt_idx]
        for(pos in seq_len(mapped)) {
          if(stats::runif(1) < mutation_prob) {
            # Current value information
            current <- ind$genotype[[nt_idx]][[pos]]
            current_depth <- current[3]

            # Generate new codon with Gaussian mutation
            new_codon <- stats::pnorm(stats::rnorm(1, stats::qnorm(current[2]), 0.5))
            new_codon <- min(max(new_codon, 0), 1)

            # Handle depth constraints
            if(current_depth >= private$.grammar$max_depth) {
              # Use shortest path rules
              shortest_path <- private$.grammar$get_non_recursive_rules(non_terminals[nt_idx])
              if(length(shortest_path$rules) > 0) {
                # Select rule based on normalized probabilities
                total_prob <- shortest_path$total_probability
                cum_prob <- 0
                for(rule in shortest_path$rules) {
                  cum_prob <- cum_prob + rule$probability/total_prob
                  if(new_codon <= cum_prob) {
                    ind$genotype[[nt_idx]][[pos]] <- c(rule$index - 1, new_codon, current_depth)
                    break
                  }
                }
              }
            } else {
              # Normal expansion using PCFG
              cum_prob <- 0
              for(rule_idx in seq_along(pcfg[nt_idx,])) {
                cum_prob <- cum_prob + pcfg[nt_idx, rule_idx]
                if(new_codon <= cum_prob) {
                  ind$genotype[[nt_idx]][[pos]] <- c(rule_idx - 1, new_codon, current_depth)
                  break
                }
              }
            }
          }
        }
      }

      ind
    },

    # Crossover operator
    .crossover = function(parent1, parent2) {
      # Create mask for selecting lists from parents
      mask <- stats::rbinom(length(parent1$genotype), 1, 0.5) == 1

      # Initialize offspring
      offspring <- list(
        genotype = vector("list", length(parent1$genotype)),
        fitness = NULL,
        mapping_values = numeric(length(parent1$genotype))
      )

      # Inherit genetic material based on mask - with deep copies
      for(i in seq_along(mask)) {
        # Deep copy the selected parent's genotype component
        parent_component <- if(mask[i]) parent1$genotype[[i]] else parent2$genotype[[i]]
        offspring$genotype[[i]] <- lapply(parent_component, function(x) c(x))  # Force copy of numeric vectors
      }

      # Compute tree depth for offspring
      phenotype_info <- private$.grammar$mapping(offspring$genotype, offspring$mapping_values)
      offspring$tree_depth <- phenotype_info$max_depth

      offspring
    },

    # Evaluate individual
    .evaluate = function(individual) {
      # Track current positions
      positions_to_map <- numeric(length(individual$genotype))

      # Map genotype to phenotype
      phenotype_info <- private$.grammar$mapping(
        mapping_rules = individual$genotype,
        positions_to_map = positions_to_map
      )

      # Update individual with mapping info
      individual$mapping_values <- positions_to_map
      individual$tree_depth <- phenotype_info$max_depth
      individual$fitness <- private$.evaluation_function(phenotype_info$phenotype)

      individual
    }
  ),

  public = list(
    #' @description
    #' Create a new PSGE instance
    #'
    #' @param grammar A Grammar object defining the solution space
    #' @param eval_func Function that takes a phenotype and returns fitness score
    #' @param population_size Integer size of population (default: 100)
    #' @param generations Integer number of generations to evolve (default: 50)
    #' @param elite_size Integer number of elite individuals (default: 10)
    #' @param tournament_size Integer size of tournament selection (default: 3)
    #' @param crossover_prob Numeric probability of crossover (default: 0.9)
    #' @param mutation_prob Numeric probability of mutation (default: 0.1)
    #' @param learning_factor Numeric PCFG learning rate (default: 0.01)
    #'
    #' @return A new `PSGE` object
    initialize = function(grammar,
                          eval_func,
                          population_size = 100,
                          generations = 50,
                          elite_size = 10,
                          tournament_size = 3,
                          crossover_prob = 0.9,
                          mutation_prob = 0.1,
                          learning_factor = 0.01) {

      private$.grammar <- grammar
      private$.evaluation_function <- eval_func
      private$.population_size <- population_size
      private$.generations <- generations
      private$.elite_size <- elite_size
      private$.tournament_size <- tournament_size
      private$.crossover_prob <- crossover_prob
      private$.mutation_prob <- mutation_prob
      private$.learning_factor <- learning_factor

      # Initialize population
      private$.population <- replicate(
        population_size,
        private$.generate_random_individual(),
        simplify = FALSE
      )
    },

    #' @description
    #' Run the evolutionary process
    #'
    #' @param verbose Logical; if TRUE, print progress information
    #' @return Best individual found (list with components genotype, fitness, etc.)
    evolve = function(verbose = TRUE) {
      # Initialize progress tracking
      generation <- 0
      update_flag <- FALSE  # Alternating flag for PCFG updates

      # Main evolution loop
      while(generation <= private$.generations) {
        # Evaluate population
        private$.population <- lapply(private$.population, private$.evaluate)

        # Sort population by fitness
        fitness_values <- sapply(private$.population, function(x) x$fitness)
        sorted_indices <- order(fitness_values)
        private$.population <- private$.population[sorted_indices]

        # Update best individual
        if(private$.population[[1]]$fitness < private$.best_fitness) {
          private$.best_individual <- utils::modifyList(list(), private$.population[[1]])
          private$.best_fitness <- private$.best_individual$fitness
        }

        # Update PCFG probabilities
        if(!update_flag) {
          # Use best overall
          private$.grammar$update_pcfg_probabilities(
            private$.best_individual,
            private$.learning_factor
          )
        } else {
          # Use best from current generation
          private$.grammar$update_pcfg_probabilities(
            private$.population[[1]],
            private$.learning_factor
          )
        }
        update_flag <- !update_flag

        # Progress report
        if(verbose && generation %% 10 == 0) {
          cat(sprintf(
            "Generation %d: Best Fitness = %.6f, Mean Fitness = %.6f\n",
            generation,
            private$.population[[1]]$fitness,
            mean(fitness_values)
          ))
        }

        # Create new population
        new_population <- vector("list", private$.population_size)

        # Elitism - remap and reevaluate elite individuals
        elite_individuals <- lapply(
          private$.population[1:private$.elite_size],
          function(ind) {
            elite <- utils::modifyList(list(), ind)  # Deep copy
            private$.evaluate(elite)  # Remap and reevaluate
          }
        )
        new_population[1:private$.elite_size] <- elite_individuals

        # Fill rest of population
        for(i in (private$.elite_size + 1):private$.population_size) {
          if(stats::runif(1) < private$.crossover_prob) {
            # Crossover
            parent1 <- private$.tournament_select()
            parent2 <- private$.tournament_select()
            offspring <- private$.crossover(parent1, parent2)
          } else {
            # Mutation only
            offspring <- private$.tournament_select()
          }

          # Apply mutation
          offspring <- private$.mutate(offspring, private$.mutation_prob)

          new_population[[i]] <- offspring
        }

        # Update population and increment generation
        private$.population <- new_population
        generation <- generation + 1
      }

      # Return best solution found
      private$.best_individual
    },

    #' @description
    #' Get the current best individual
    #'
    #' @return List containing the best individual's information
    get_best = function() {
      private$.best_individual
    },

    #' @description
    #' Get the current population
    #'
    #' @return List of individuals in the current population
    get_population = function() {
      private$.population
    }
  )
)

#' Create a new PSGE instance
#'
#' @description
#' Creates a new Probabilistic Structured Grammatical Evolution (PSGE) instance
#' with the specified parameters. This is the recommended way to create a new
#' PSGE object rather than calling the R6 constructor directly.
#'
#' @param grammar A Grammar object defining the solution space
#' @param eval_func Function that takes a phenotype string and returns a numeric
#'   fitness value to be minimized
#' @param population_size Integer size of population (default: 100)
#' @param generations Integer number of generations to evolve (default: 50)
#' @param elite_size Integer number of elite individuals to preserve (default: 10)
#' @param tournament_size Integer size of tournament selection (default: 3)
#' @param crossover_prob Numeric probability of crossover [0,1] (default: 0.9)
#' @param mutation_prob Numeric probability of mutation [0,1] (default: 0.1)
#' @param learning_factor Numeric PCFG learning rate (0,1] (default: 0.01)
#'
#' @return A new PSGE object
#'
#' @examples
#' # Create a simple grammar
#' grammar_text <- "
#' <expr> ::= <var> | <expr> <op> <expr>
#' <op> ::= + | -
#' <var> ::= x | 1.0"
#' grammar_file <- tempfile(fileext = ".bnf")
#' writeLines(grammar_text, grammar_file)
#' grammar <- create_grammar(grammar_file)
#'
#' # Define a simple evaluation function
#' eval_func <- function(phenotype) {
#'   # Return mean squared error from f(x) = x + 1
#'   x <- seq(-1, 1, by=0.1)
#'   y_true <- x + 1
#'   y_pred <- eval(parse(text=phenotype))
#'   mean((y_true - y_pred)^2)
#' }
#'
#' # Create PSGE instance
#' psge <- create_psge(
#'   grammar = grammar,
#'   eval_func = eval_func,
#'   population_size = 50,
#'   generations = 20
#' )
#'
#' @export
create_psge <- function(grammar,
                        eval_func,
                        population_size = 100,
                        generations = 50,
                        elite_size = 10,
                        tournament_size = 3,
                        crossover_prob = 0.9,
                        mutation_prob = 0.1,
                        learning_factor = 0.01) {

  # Validate inputs
  if (!inherits(grammar, "Grammar")) {
    stop("'grammar' must be a Grammar object")
  }
  if (!is.function(eval_func)) {
    stop("'eval_func' must be a function")
  }
  if (!is.numeric(population_size) || population_size < 1) {
    stop("'population_size' must be a positive integer")
  }
  if (!is.numeric(generations) || generations < 1) {
    stop("'generations' must be a positive integer")
  }
  if (!is.numeric(elite_size) || elite_size < 0 || elite_size >= population_size) {
    stop("'elite_size' must be a non-negative integer less than population_size")
  }
  if (!is.numeric(tournament_size) || tournament_size < 1 ||
      tournament_size > population_size) {
    stop("'tournament_size' must be an integer between 1 and population_size")
  }
  if (!is.numeric(crossover_prob) || crossover_prob < 0 || crossover_prob > 1) {
    stop("'crossover_prob' must be between 0 and 1")
  }
  if (!is.numeric(mutation_prob) || mutation_prob < 0 || mutation_prob > 1) {
    stop("'mutation_prob' must be between 0 and 1")
  }
  if (!is.numeric(learning_factor) || learning_factor <= 0 || learning_factor > 1) {
    stop("'learning_factor' must be in (0,1]")
  }

  # Create PSGE instance with validated parameters
  PSGE$new(
    grammar = grammar,
    eval_func = eval_func,
    population_size = as.integer(population_size),
    generations = as.integer(generations),
    elite_size = as.integer(elite_size),
    tournament_size = as.integer(tournament_size),
    crossover_prob = crossover_prob,
    mutation_prob = mutation_prob,
    learning_factor = learning_factor
  )
}

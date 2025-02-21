#' Run the PSGE algorithm with progress reporting
#'
#' @param grammar PCFG grammar structure from read_pcfg_grammar()
#' @param fitness_fn Fitness evaluation function
#' @param pop_size Population size
#' @param generations Number of generations
#' @param crossover_rate Probability of crossover
#' @param mutation_rate Mutation rate
#' @param max_depth Maximum tree depth
#' @param tournament_size Tournament selection size
#' @param elitism Number of elite individuals to preserve
#' @param learning_factor Learning rate for grammar probability updates
#' @param alternate_update Whether to alternate between best overall and best in generation
#' @param force_perturbation Whether to force perturbation of all probability distributions
#' @param perturbation_factor Amount of perturbation when forcing changes
#' @param verbose Whether to print additional progress information
#' @param parallel Whether to use parallel evaluation
#' @param n_workers Number of workers for parallel processing
#' @param future_strategy Future strategy to use for parallel processing
#' @param post_generation_fn Function to run after population generation (optional)
#' @param post_evaluation_fn Function to run after fitness evaluation (optional)
#'
#' @details
#' Note that we use progressr progress bars with this function:
#'
#' ```r
#' library(progressr)
#'
#' # Configure progress reporting
#' handlers(global = TRUE)
#' handlers("progress")
#'
#' # Run with progress reporting
#' with_progress({
#'   result <- run_psge(
#'     grammar = my_grammar,
#'     fitness_fn = my_fitness_function,
#'     generations = 100
#'   )
#' })
#' ```
#'
#' @return List containing best solution and final grammar
#' @export
run_psge <- function(grammar, fitness_fn, pop_size = 100, generations = 50,
                     crossover_rate = 0.9, mutation_rate = 0.1, max_depth = 17,
                     tournament_size = 3, elitism = NULL, learning_factor = 0.01,
                     alternate_update = TRUE, force_perturbation = FALSE,
                     perturbation_factor = 0.1, verbose = TRUE,
                     parallel = FALSE, n_workers = 2, future_strategy = "multisession",
                     post_generation_fn = NULL, post_evaluation_fn = NULL) {
  # Create progress handlers
  p <- progressr::progressor(steps = generations + 2) # +2 for initialization and final processing

  # Set default elitism if not specified
  if (is.null(elitism)) {
    elitism <- max(round(pop_size * 0.1), 1)
  }

  # Create a deep copy of the grammar to avoid modifying the original
  working_grammar <- grammar

  # Initialize population
  if (verbose) p("Initializing population...")
  population <- initialize_population(working_grammar, pop_size, max_depth)

  # Create initial context
  context <- list(
    generation = 0,
    total_generations = generations,
    population_size = pop_size,
    current_best_fitness = Inf,
    best_ever_fitness = Inf,
    mean_fitness = NA,
    elapsed_time = 0,
    start_time = Sys.time(),
    grammar = working_grammar
  )

  # Run post-generation function if provided
  if (!is.null(post_generation_fn)) {
    population <- post_generation_fn(population, context)
    if (length(population) != pop_size) {
      warning("post_generation_fn changed population size. Adjusting back to pop_size.")
      if (length(population) > pop_size) {
        population <- population[1:pop_size]
      } else {
        additional <- initialize_population(
          working_grammar,
          pop_size - length(population),
          max_depth
        )
        population <- c(population, additional)
      }
    }
  }

  # Evaluate initial population
  if (verbose) p("Evaluating initial population...")
  population <- evaluate_population(
    population, working_grammar, fitness_fn, max_depth,
    parallel = parallel, n_workers = n_workers,
    future_strategy = future_strategy
  )

  # Update context with initial fitness information
  fitness_values <- sapply(population, function(ind) ind$fitness)
  context$current_best_fitness <- min(fitness_values)
  context$best_ever_fitness <- context$current_best_fitness
  context$mean_fitness <- mean(fitness_values)
  context$elapsed_time <- as.numeric(difftime(Sys.time(), context$start_time, units = "secs"))

  # Run post-evaluation function if provided
  if (!is.null(post_evaluation_fn)) {
    population <- post_evaluation_fn(population, context)
    if (length(population) != pop_size) {
      warning("post_evaluation_fn changed population size. Adjusting back to pop_size.")
      if (length(population) > pop_size) {
        population <- population[1:pop_size]
      } else {
        additional <- initialize_population(
          working_grammar,
          pop_size - length(population),
          max_depth
        )
        additional <- evaluate_population(
          additional, working_grammar, fitness_fn, max_depth,
          parallel = parallel, n_workers = n_workers,
          future_strategy = future_strategy
        )
        population <- c(population, additional)
      }
    }
  }

  # Sort population
  population <- population[order(sapply(population, function(ind) ind$fitness))]

  # Initialize tracking variables
  best_overall <- population[[1]]
  history <- data.frame(
    generation = 0,
    best_fitness = best_overall$fitness,
    mean_fitness = mean(sapply(population, function(ind) ind$fitness))
  )

  # Main evolutionary loop
  for (gen in 1:generations) {
    # Update progress with generation info
    progress_msg <- sprintf("Generation %d/%d - Best fitness: %.6f",
                            gen, generations, population[[1]]$fitness)
    p(progress_msg)

    # Update context for this generation
    context$generation <- gen
    context$current_best_fitness <- population[[1]]$fitness
    context$best_ever_fitness <- min(context$best_ever_fitness, context$current_best_fitness)
    context$mean_fitness <- mean(sapply(population, function(ind) ind$fitness))
    context$elapsed_time <- as.numeric(difftime(Sys.time(), context$start_time, units = "secs"))
    context$grammar <- working_grammar

    # Update best individuals
    current_best <- population[[1]]
    if (current_best$fitness < best_overall$fitness) {
      best_overall <- current_best
    }

    # Update grammar probabilities
    best_to_update <- if (alternate_update && gen %% 2 == 0) best_overall else current_best
    working_grammar <- update_grammar_probabilities(
      working_grammar, best_to_update, learning_factor,
      force_perturbation, perturbation_factor
    )

    # Create new generation
    population <- create_next_generation(
      population, working_grammar, crossover_rate, mutation_rate,
      tournament_size, elitism
    )

    # Run post-generation function if provided
    if (!is.null(post_generation_fn)) {
      population <- post_generation_fn(population, context)
      if (length(population) != pop_size) {
        warning("post_generation_fn changed population size in generation ", gen, ". Adjusting.")
        if (length(population) > pop_size) {
          population <- population[1:pop_size]
        } else {
          additional <- initialize_population(working_grammar, pop_size - length(population), max_depth)
          population <- c(population, additional)
        }
      }
    }

    # Evaluate new population
    population <- evaluate_population(
      population, working_grammar, fitness_fn, max_depth,
      parallel = parallel, n_workers = n_workers,
      future_strategy = future_strategy
    )

    # Run post-evaluation function if provided
    if (!is.null(post_evaluation_fn)) {
      population <- post_evaluation_fn(population, context)
      if (length(population) != pop_size) {
        warning("post_evaluation_fn changed population size in generation ", gen, ". Adjusting.")
        if (length(population) > pop_size) {
          population <- population[1:pop_size]
        } else {
          additional <- initialize_population(working_grammar, pop_size - length(population), max_depth)
          additional <- evaluate_population(
            additional, working_grammar, fitness_fn, max_depth,
            parallel = parallel, n_workers = n_workers,
            future_strategy = future_strategy
          )
          population <- c(population, additional)
        }
      }
    }

    # Sort population
    population <- population[order(sapply(population, function(ind) ind$fitness))]

    # Track history
    history <- rbind(history, data.frame(
      generation = gen,
      best_fitness = population[[1]]$fitness,
      mean_fitness = mean(sapply(population, function(ind) ind$fitness))
    ))
  }

  # Final progress update
  p("Evolution completed")

  if (verbose) {
    cat(sprintf("\nBest fitness: %f\n", best_overall$fitness))
    cat("Best phenotype:", best_overall$phenotype, "\n")
  }

  return(list(
    best_solution = best_overall,
    final_grammar = working_grammar,
    history = history
  ))
}

#' Complete PSGE pipeline
#'
#' @param grammar_file Path to BNF grammar file
#' @param fitness_fn Fitness evaluation function
#' @param ... Additional parameters passed to run_psge
#' @return Best solution and final grammar
#' @export
psge <- function(grammar_file, fitness_fn, ...) {
  # Load and parse grammar
  grammar <- read_pcfg_grammar(grammar_file)

  # Run PSGE algorithm
  result <- run_psge(grammar, fitness_fn, ...)

  return(result)
}

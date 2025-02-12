#' Create and Configure PSGE Instance
#'
#' @param grammar_path Path to BNF grammar file
#' @param fitness_fn Function that evaluates phenotypes
#' @param pop_size Population size (default: 100)
#' @param generations Number of generations (default: 50)
#' @param crossover_rate Probability of crossover (default: 0.9)
#' @param mutation_rate Probability of mutation (default: 0.1)
#' @param max_depth Maximum tree depth (default: 17)
#' @param learning_rate PCFG learning rate (default: 0.01)
#' @param tournament_size Tournament selection size (default: 3)
#' @param elitism Number of elites to preserve (default: NULL, 10% of pop_size)
#' @param adaptive_learning Enables adaptive learning factor (default: FALSE)
#' @param seed Random seed for reproducibility (default: NULL)
#' @return Configured PSGE_EA instance
#' @export
setup_psge <- function(grammar_path,
                       fitness_fn,
                       pop_size = 100,
                       generations = 50,
                       crossover_rate = 0.9,
                       mutation_rate = 0.1,
                       max_depth = 17,
                       learning_rate = 0.01,
                       tournament_size = 3,
                       elitism = NULL,
                       adaptive_learning = FALSE,
                       seed = NULL) {

  # Input validation
  if (!file.exists(grammar_path)) {
    stop("Grammar file not found: ", grammar_path)
  }
  if (!is.function(fitness_fn)) {
    stop("fitness_fn must be a function")
  }

  # Create and initialize grammar
  grammar <- Grammar$new()
  grammar$set_path(grammar_path)
  grammar$set_max_tree_depth(max_depth)
  grammar$read_grammar()

  # Set default elitism if not specified (10% of population)
  if (is.null(elitism)) {
    elitism <- max(round(pop_size * 0.1), 1)
  }

  # Configure parameters
  params <- list(
    popsize = pop_size,
    generations = generations,
    prob_crossover = crossover_rate,
    prob_mutation = mutation_rate,
    learning_factor = learning_rate,
    tournament_size = tournament_size,
    elitism = elitism,
    adaptive_lf = adaptive_learning,
    seed = seed
  )

  # Create EA instance
  ea <- Evolution$new(
    grammar = grammar,
    eval_func = fitness_fn,
    params = params
  )

  return(ea)
}

#' Run PSGE Evolution
#'
#' @param psge_instance PSGE instance from setup_psge()
#' @param verbose Print progress (default TRUE)
#' @param progress_fn Custom progress callback function (optional)
#' @return Best solution found
#' @export
run_psge <- function(psge_instance,
                     verbose = TRUE,
                     progress_fn = NULL) {

  if (!inherits(psge_instance, "Evolution")) {
    stop("psge_instance must be a Evolution object")
  }

  # Default progress function
  default_progress <- function(gen, fitness) {
    cat(sprintf("Generation %d: Best Fitness = %f\n", gen, fitness))
  }

  # Set progress callback
  if (verbose && is.null(progress_fn)) {
    progress_fn <- default_progress
  }

  # Store original output settings
  original_messages <- getOption("warning.length")
  options(warning.length = 8170)

  # Run evolution with progress tracking
  tryCatch({
    solution <- psge_instance$run_evolution()

    if (verbose) {
      cat("\nEvolution completed\n")
      cat(sprintf("Best fitness: %f\n", solution$fitness))
      cat("Best phenotype:", solution$phenotype, "\n")
    }

    return(solution)

  }, error = function(e) {
    stop("Evolution failed: ", e$message)
  }, finally = {
    # Restore original settings
    options(warning.length = original_messages)
  })
}

#' Complete PSGE Pipeline
#'
#' @param grammar_path Path to BNF grammar file
#' @param fitness_fn Fitness function
#' @param ... Additional parameters passed to setup_psge()
#' @return Best solution found
#' @export
psge <- function(grammar_path, fitness_fn, ...) {
  instance <- setup_psge(grammar_path, fitness_fn, ...)
  return(run_psge(instance))
}

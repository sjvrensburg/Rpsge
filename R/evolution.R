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

#' Initialize a Population with Both Directed and Random Individuals Using Enhanced Parallel Processing
#'
#' Creates a population of size \code{pop_size} using parallel processing with improved
#' duplicate handling across concurrent tasks. A portion of individuals is generated
#' using directed templates from GADI, while the rest are generated randomly.
#'
#' @param grammar A PCFG grammar structure (e.g., from \code{read_pcfg_grammar()}).
#' @param pop_size Desired population size (integer).
#' @param max_depth Max derivation depth to prevent overly large expansions.
#' @param use_gadi Logical indicating whether to use directed genotype creation
#'   (TRUE) or just random initialization (FALSE).
#' @param fraction_directed Fraction (0 to 1) of individuals to create via
#'   directed templates. Default 0.3 means 30% of the population is
#'   template-based, and 70% is random.
#' @param analysis If already analyzed grammar diversity is available, pass it
#'   here. Otherwise, if \code{use_gadi} is TRUE and no analysis is provided, we
#'   will call \code{analyze_grammar_for_diversity()} internally.
#' @param check_duplicates Logical indicating whether to avoid duplicate phenotypes
#'   in the initial population. Default TRUE.
#' @param max_attempts Maximum number of attempts to generate a unique individual.
#'   Used only if check_duplicates is TRUE. Default 100.
#' @param parallel Logical: whether to use parallel processing. Default TRUE.
#' @param n_cores Integer: number of CPU cores to use for parallel processing.
#'   If NULL (default), it uses all available cores minus 1.
#' @param batch_size Integer: number of individuals to generate in each parallel batch.
#'   Default 10. Smaller values reduce collision risk but increase overhead.
#' @param sequential_batches Logical: whether to process batches sequentially but with
#'   internal parallelism (TRUE) or fully in parallel (FALSE). Default TRUE for better
#'   duplicate control.
#' @param post_process_duplicates Logical: whether to filter out any duplicates in a
#'   final post-processing step. Default TRUE.
#' @param exact_population_size Logical: whether to ensure exactly pop_size individuals
#'   are returned (TRUE) or accept slightly fewer if max_attempts is reached (FALSE).
#'   Default TRUE.
#' @param verbose Logical indicating whether to print information about the
#'   initialization process. Default FALSE.
#' @param ... Additional arguments passed to \code{directed_templates_from_analysis()}.
#'
#' @return A list of individuals, where each individual is a list containing:
#'   \itemize{
#'     \item genotype: Named list of numeric vectors (codons for each non-terminal).
#'     \item phenotype: The derived string or expression from mapping the genotype.
#'     \item fitness: Initialized to NULL here (not yet evaluated).
#'     \item mapping_positions: Additional info from \code{map_genotype()}.
#'     \item tree_depth: Depth of the derived parse tree.
#'     \item initialization: Character indicating how the individual was created ("directed" or "random").
#'   }
#'
#' @examples
#' \dontrun{
#'   # Basic usage with parallel processing
#'   library(future)
#'   library(future.apply)
#'
#'   grammar <- read_pcfg_grammar("path/to/grammar.bnf")
#'   pop <- initialize_population(
#'     grammar, 100,
#'     use_gadi = TRUE,
#'     fraction_directed = 0.3,
#'     parallel = TRUE,
#'     sequential_batches = TRUE  # Process batches sequentially with internal parallelism
#'   )
#'
#'   # Maximum parallelism but with potential duplicate risk
#'   pop2 <- initialize_population(
#'     grammar, 200,
#'     parallel = TRUE,
#'     sequential_batches = FALSE,  # Fully parallel processing
#'     post_process_duplicates = TRUE  # Remove any duplicates at the end
#'   )
#' }
#'
#' @export
initialize_population <- function(
    grammar,
    pop_size,
    max_depth = 17,
    use_gadi = TRUE,
    fraction_directed = 0.3,
    analysis = NULL,
    check_duplicates = TRUE,
    max_attempts = 100,
    parallel = TRUE,
    n_cores = NULL,
    batch_size = 10,
    sequential_batches = TRUE,
    post_process_duplicates = TRUE,
    exact_population_size = TRUE,
    verbose = FALSE,
    ...
) {
  # Validate input parameters
  if (pop_size < 1) stop("Population size must be at least 1")
  if (fraction_directed < 0 || fraction_directed > 1) {
    stop("fraction_directed must be between 0 and 1")
  }

  # Check for required packages
  if (parallel && !requireNamespace("future", quietly = TRUE)) {
    warning("The 'future' package is required for parallel processing. Falling back to sequential mode.")
    parallel <- FALSE
  }
  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    warning("The 'future.apply' package is required for parallel processing. Falling back to sequential mode.")
    parallel <- FALSE
  }

  start_time <- Sys.time()

  # Set up parallel backend if requested
  if (parallel) {
    # Store the current plan to restore it later
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    # Determine number of cores
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }

    if (verbose) {
      cat("Setting up parallel processing with", n_cores, "cores\n")
    }

    # Set up multisession parallel backend
    future::plan(future::multisession, workers = n_cores)
  }

  # Create shared environment for concurrent duplicate tracking if needed
  create_shared_phenotype_tracker <- function() {
    shared_env <- new.env()
    shared_env$phenotypes <- character(0)
    shared_env$lock_count <- 0  # For diagnostic purposes
    shared_env$locked <- FALSE

    # Create lock acquisition function with timeout
    shared_env$acquire_lock <- function(timeout = 5) {
      start_time <- Sys.time()
      while (shared_env$locked) {
        # Check for timeout to prevent deadlocks
        if (difftime(Sys.time(), start_time, units = "secs") > timeout) {
          stop("Timeout waiting for lock acquisition")
        }
        Sys.sleep(0.001)  # Small sleep to reduce CPU usage
      }
      shared_env$locked <- TRUE
      shared_env$lock_count <- shared_env$lock_count + 1
      return(TRUE)
    }

    # Create lock release function
    shared_env$release_lock <- function() {
      shared_env$locked <- FALSE
      return(TRUE)
    }

    # Safe operation that automatically handles locking and unlocking
    shared_env$with_lock <- function(expr) {
      locked <- FALSE
      tryCatch({
        # Acquire lock
        shared_env$acquire_lock()
        locked <- TRUE

        # Evaluate expression in the locked context
        result <- eval(expr, envir = parent.frame())

        # Release lock
        shared_env$release_lock()
        locked <- FALSE

        return(result)
      },
      error = function(e) {
        # Make sure we release the lock if an error occurs
        if (locked) shared_env$release_lock()
        stop(e)
      },
      finally = {
        # Extra safety to ensure lock is released even in exceptional cases
        if (locked) shared_env$release_lock()
      })
    }

    # Function to safely check and add a phenotype
    # Returns TRUE if the phenotype was added (not a duplicate)
    # Returns FALSE if the phenotype already exists
    shared_env$check_and_add_phenotype <- function(phenotype) {
      # Use the with_lock pattern to ensure safety
      shared_env$with_lock({
        if (phenotype %in% shared_env$phenotypes) {
          return(FALSE)  # Duplicate found
        } else {
          # Add to shared list
          shared_env$phenotypes <- c(shared_env$phenotypes, phenotype)
          return(TRUE)   # New phenotype added
        }
      })
    }

    # Function to get current unique phenotype count
    shared_env$get_phenotype_count <- function() {
      shared_env$with_lock({
        return(length(shared_env$phenotypes))
      })
    }

    # Function to check if a phenotype exists without adding it
    shared_env$check_phenotype <- function(phenotype) {
      shared_env$with_lock({
        return(phenotype %in% shared_env$phenotypes)
      })
    }

    return(shared_env)
  }

  # Set up shared environment for concurrent duplicate tracking if needed
  shared_env <- NULL
  if (parallel && !sequential_batches && check_duplicates) {
    shared_env <- create_shared_phenotype_tracker()

    if (verbose) {
      cat("Using concurrent phenotype tracking with safe locking pattern\n")
    }
  }

  # If user wants GADI but didn't provide an analysis, do it now
  if (use_gadi && is.null(analysis)) {
    if (verbose) cat("Analyzing grammar diversity...\n")
    # Extract GADI-specific arguments
    gadi_args <- list(...)
    analysis_args <- list(grammar = grammar)

    # Only include parameters that analyze_grammar_for_diversity can accept
    valid_analysis_params <- c("use_semantic_awareness", "semantic_group_fn")
    for (param in intersect(names(gadi_args), valid_analysis_params)) {
      analysis_args[[param]] <- gadi_args[[param]]
    }

    # Perform the analysis with appropriate arguments
    analysis <- do.call(analyze_grammar_for_diversity, analysis_args)
  }

  # Generate directed templates if GADI is on
  templates <- list()
  if (use_gadi && !is.null(analysis)) {
    if (verbose) cat("Generating directed templates...\n")
    templates <- directed_templates_from_analysis(analysis, ...)

    if (length(templates) == 0 && verbose) {
      warning("No templates were generated. Check your grammar analysis or parameters.")
    } else if (verbose) {
      cat("Generated", length(templates), "directed templates\n")
    }
  }

  # Figure out how many individuals should come from templates
  directed_count <- round(pop_size * fraction_directed)
  random_count <- pop_size - directed_count

  if (verbose) {
    cat("Creating population of size", pop_size, "\n")
    cat("  -", directed_count, "directed individuals\n")
    cat("  -", random_count, "random individuals\n")
  }

  # Helper function to check for duplicate phenotypes
  is_duplicate_phenotype <- function(new_pheno, existing_phenos) {
    if (length(existing_phenos) == 0) return(FALSE)
    return(new_pheno %in% existing_phenos)
  }

  # Helper function to generate a single individual
  generate_individual <- function(method, template = NULL, existing_phenos = NULL, max_attempts = 100) {
    for (attempt in 1:max_attempts) {
      if (method == "directed" && !is.null(template)) {
        # Create individual from template
        indiv <- create_individual_from_template(grammar, template, max_depth = max_depth)
        indiv$initialization <- "directed"
      } else {
        # Create random individual
        rand_genotype <- generate_random_individual(grammar, max_depth)
        mapping_result <- map_genotype(grammar, rand_genotype, max_depth)

        indiv <- list(
          genotype = rand_genotype,
          phenotype = mapping_result$phenotype,
          fitness = NULL,
          mapping_positions = mapping_result$positions,
          tree_depth = mapping_result$depth,
          initialization = "random"
        )
      }

      # Check for duplicate phenotype if required
      if (!check_duplicates || !is_duplicate_phenotype(indiv$phenotype, existing_phenos)) {
        return(indiv)  # Found a non-duplicate
      }
    }

    # If we reach here, max_attempts were exceeded
    return(NULL)
  }

  # Function to generate a batch of individuals
  generate_batch <- function(batch_spec, existing_phenotypes = NULL) {
    results <- list()

    # If we're in parallel mode and using the shared environment for concurrent tracking
    if (!is.null(shared_env) && is.environment(existing_phenotypes) &&
        exists("check_and_add_phenotype", envir = existing_phenotypes)) {

      # Use the shared environment for concurrent duplicate tracking
      for (spec in batch_spec) {
        # Try to generate a unique individual
        indiv <- NULL
        attempt <- 0

        while (is.null(indiv) && attempt < max_attempts) {
          attempt <- attempt + 1

          # Generate candidate
          if (spec$method == "directed" && !is.null(spec$template)) {
            candidate <- create_individual_from_template(grammar, spec$template, max_depth = max_depth)
            candidate$initialization <- "directed"
          } else {
            rand_genotype <- generate_random_individual(grammar, max_depth)
            mapping_result <- map_genotype(grammar, rand_genotype, max_depth)

            candidate <- list(
              genotype = rand_genotype,
              phenotype = mapping_result$phenotype,
              fitness = NULL,
              mapping_positions = mapping_result$positions,
              tree_depth = mapping_result$depth,
              initialization = "random"
            )
          }

          # Use the encapsulated function for thread-safe checking and adding
          if (existing_phenotypes$check_and_add_phenotype(candidate$phenotype)) {
            indiv <- candidate  # Not a duplicate, keep it
          }
        }

        if (!is.null(indiv)) {
          results <- c(results, list(indiv))
        }
      }
    } else {
      # Standard mode - check against local list
      current_phenotypes <- existing_phenotypes

      for (spec in batch_spec) {
        indiv <- generate_individual(
          method = spec$method,
          template = spec$template,
          existing_phenos = current_phenotypes,
          max_attempts = max_attempts
        )

        if (!is.null(indiv)) {
          results <- c(results, list(indiv))

          # Add this phenotype to the list of existing ones
          if (check_duplicates) {
            current_phenotypes <- c(current_phenotypes, indiv$phenotype)
          }
        }
      }
    }

    return(results)
  }

  # Prepare template specifications for directed individuals
  directed_specs <- list()

  if (directed_count > 0 && length(templates) > 0) {
    # If we have more templates than needed, randomly sample
    if (length(templates) > directed_count) {
      template_indices <- sample(length(templates), directed_count)
      used_templates <- templates[template_indices]
    } else {
      # If we need more directed individuals than templates,
      # we'll reuse templates (cycling through them)
      cycles_needed <- ceiling(directed_count / length(templates))
      template_indices <- rep(1:length(templates), times = cycles_needed)
      template_indices <- template_indices[1:directed_count]
      used_templates <- templates[template_indices]
    }

    directed_specs <- lapply(used_templates, function(tmpl) {
      list(method = "directed", template = tmpl)
    })
  }

  # Prepare random specifications
  random_specs <- replicate(random_count, list(method = "random", template = NULL), simplify = FALSE)

  # Combine all specifications
  all_specs <- c(directed_specs, random_specs)

  # Final population
  population <- list()
  existing_phenotypes <- character(0)

  # Split into batches for processing
  n_batches <- ceiling(length(all_specs) / batch_size)
  batch_specs <- split(all_specs, ceiling(seq_along(all_specs) / batch_size))

  if (verbose) {
    cat("Processing", length(batch_specs), "batches with batch size", batch_size, "\n")
  }

  if (parallel) {
    if (sequential_batches) {
      # Process batches sequentially, but with internal parallelism
      if (verbose) cat("Processing batches sequentially with internal parallelism...\n")

      for (batch_idx in seq_along(batch_specs)) {
        batch <- batch_specs[[batch_idx]]

        if (verbose) {
          cat("  Processing batch", batch_idx, "of", length(batch_specs),
              "(", length(batch), "individuals )...\n")
        }

        # Process this batch in parallel
        batch_results <- future.apply::future_lapply(
          batch,
          function(spec) {
            generate_individual(
              method = spec$method,
              template = spec$template,
              existing_phenos = existing_phenotypes,
              max_attempts = max_attempts
            )
          },
          future.seed = TRUE
        )

        # Remove NULLs (failed attempts)
        batch_results <- batch_results[!sapply(batch_results, is.null)]

        # Add to population
        population <- c(population, batch_results)

        # Update existing phenotypes for the next batch
        if (check_duplicates) {
          existing_phenotypes <- c(
            existing_phenotypes,
            sapply(batch_results, function(ind) ind$phenotype)
          )
        }

        if (verbose) {
          cat("    Added", length(batch_results), "individuals from batch", batch_idx, "\n")
          cat("    Current population size:", length(population), "\n")
        }

        # Stop if we've reached the desired population size
        if (length(population) >= pop_size) {
          if (verbose) cat("  Target population size reached. Stopping batch processing.\n")
          break
        }
      }
    } else {
      # Process all batches fully in parallel
      if (verbose) cat("Processing all batches in parallel...\n")

      # If using shared environment for concurrent duplicate checking
      if (!is.null(shared_env)) {
        # First add any existing phenotypes to the shared environment
        if (length(existing_phenotypes) > 0) {
          for (pheno in existing_phenotypes) {
            shared_env$check_and_add_phenotype(pheno)
          }
        }

        # Process batches with shared duplicate tracking
        results <- future.apply::future_lapply(
          batch_specs,
          function(batch) generate_batch(batch, shared_env),
          future.seed = TRUE
        )
      } else {
        # Process batches without concurrent duplicate checking
        # (Each batch only checks against the initial phenotypes)
        pheno_list <- replicate(length(batch_specs), existing_phenotypes, simplify = FALSE)

        results <- future.apply::future_mapply(
          generate_batch,
          batch_spec = batch_specs,
          existing_phenotypes = pheno_list,
          SIMPLIFY = FALSE,
          future.seed = TRUE
        )
      }

      # Combine results into population
      for (batch_result in results) {
        population <- c(population, batch_result)
      }
    }
  } else {
    # Sequential processing (no parallelism)
    if (verbose) cat("Generating population sequentially...\n")

    # Process batches one by one
    for (batch_idx in seq_along(batch_specs)) {
      batch <- batch_specs[[batch_idx]]

      if (verbose && batch_idx %% 5 == 0) {
        cat("  Processing batch", batch_idx, "of", length(batch_specs), "\n")
      }

      batch_results <- generate_batch(batch, existing_phenotypes)

      # Add to population
      population <- c(population, batch_results)

      # Update existing phenotypes for the next batch
      if (check_duplicates) {
        existing_phenotypes <- c(
          existing_phenotypes,
          sapply(batch_results, function(ind) ind$phenotype)
        )
      }

      # Stop if we've reached the desired population size
      if (length(population) >= pop_size) {
        break
      }
    }
  }

  # Post-processing: Remove any duplicates that might have slipped through
  if (post_process_duplicates && check_duplicates) {
    if (verbose) cat("Post-processing: Checking for duplicates...\n")

    # Get unique phenotypes
    unique_phenos <- character(0)
    unique_indices <- integer(0)

    for (i in seq_along(population)) {
      pheno <- population[[i]]$phenotype
      if (!pheno %in% unique_phenos) {
        unique_phenos <- c(unique_phenos, pheno)
        unique_indices <- c(unique_indices, i)
      }
    }

    # If duplicates were found, keep only unique individuals
    if (length(unique_indices) < length(population)) {
      if (verbose) {
        cat("  Found", length(population) - length(unique_indices),
            "duplicate individuals. Removing them.\n")
      }
      population <- population[unique_indices]
    }
  }

  # Ensure exact population size if requested
  if (exact_population_size && length(population) != pop_size) {
    if (length(population) < pop_size) {
      # Need to add more individuals
      additional_needed <- pop_size - length(population)

      if (verbose) {
        cat("Need", additional_needed, "more individuals to reach target population size.\n")
      }

      # Try to add more random individuals
      additional_attempts <- 0
      max_additional_attempts <- additional_needed * max_attempts

      while (length(population) < pop_size && additional_attempts < max_additional_attempts) {
        additional_attempts <- additional_attempts + 1

        # Generate a random individual
        rand_genotype <- generate_random_individual(grammar, max_depth)
        mapping_result <- map_genotype(grammar, rand_genotype, max_depth)

        rand_ind <- list(
          genotype = rand_genotype,
          phenotype = mapping_result$phenotype,
          fitness = NULL,
          mapping_positions = mapping_result$positions,
          tree_depth = mapping_result$depth,
          initialization = "random_additional"
        )

        # Check for duplicates if requested
        if (!check_duplicates || !rand_ind$phenotype %in% sapply(population, function(ind) ind$phenotype)) {
          population <- c(population, list(rand_ind))
        }
      }

      if (length(population) < pop_size && verbose) {
        warning("Could only generate ", length(population), " unique individuals out of ",
                pop_size, " requested after additional attempts.")
      }
    } else if (length(population) > pop_size) {
      # Need to remove some individuals
      if (verbose) {
        cat("Trimming population from", length(population), "to", pop_size, "individuals.\n")
      }

      # Keep a balanced proportion of directed vs random individuals
      directed_indices <- which(sapply(population, function(ind) ind$initialization == "directed"))
      random_indices <- which(sapply(population, function(ind) ind$initialization != "directed"))

      target_directed <- min(length(directed_indices), round(pop_size * fraction_directed))
      target_random <- pop_size - target_directed

      # Random sampling to maintain diversity
      if (length(directed_indices) > target_directed) {
        directed_indices <- sample(directed_indices, target_directed)
      }

      if (length(random_indices) > target_random) {
        random_indices <- sample(random_indices, target_random)
      }

      # Combine indices and extract final population
      final_indices <- c(directed_indices, random_indices)
      population <- population[final_indices]
    }
  }

  # Provide population statistics
  if (verbose) {
    end_time <- Sys.time()
    time_taken <- difftime(end_time, start_time, units = "secs")

    cat("\nPopulation initialization complete.\n")
    cat("Final population size:", length(population), "\n")
    cat("Time taken:", round(time_taken, 2), "seconds\n")

    if (length(population) > 0) {
      # Count unique phenotypes
      unique_phenos <- unique(sapply(population, function(ind) ind$phenotype))
      cat("Unique phenotypes:", length(unique_phenos), "\n")
      cat("Effective diversity:", round(length(unique_phenos) / length(population), 2), "\n")

      # Count directed vs random
      directed_final <- sum(sapply(population, function(ind) ind$initialization == "directed"))
      random_final <- sum(sapply(population, function(ind) ind$initialization != "directed"))
      cat("Final mix:", directed_final, "directed,", random_final, "random\n")
    }
  }

  # Return the final population
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

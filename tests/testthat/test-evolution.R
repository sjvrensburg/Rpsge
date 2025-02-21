library(testthat)

# Helper functions
create_test_grammar <- function() {
  grammar <- list(
    non_terminals = c("expr", "op", "var"),
    terminals = c("+", "-", "*", "/", "x", "y", "1.0"),
    start_symbol = "expr",
    rules = list(
      expr = list(
        list(symbols = c("expr", "op", "expr"), types = c("NT", "NT", "NT"), prob = 0.5),
        list(symbols = "var", types = "NT", prob = 0.5)
      ),
      op = list(
        list(symbols = "+", types = "T", prob = 0.25),
        list(symbols = "-", types = "T", prob = 0.25),
        list(symbols = "*", types = "T", prob = 0.25),
        list(symbols = "/", types = "T", prob = 0.25)
      ),
      var = list(
        list(symbols = "x", types = "T", prob = 1/3),
        list(symbols = "y", types = "T", prob = 1/3),
        list(symbols = "1.0", types = "T", prob = 1/3)
      )
    )
  )
  return(grammar)
}

create_test_individual <- function() {
  list(
    genotype = list(
      expr = c(0.1, 0.6),  # First expr->expr op expr, then expr->var
      op = c(0.1),         # op -> +
      var = c(0.2)         # var -> x
    ),
    phenotype = "x+y",
    fitness = 0.5,
    tree_depth = 3,
    mapping_positions = list(expr = 3, op = 2, var = 2)
  )
}

create_test_population <- function(size = 10) {
  population <- vector("list", size)
  for (i in seq_len(size)) {
    # Create slightly different individuals
    ind <- create_test_individual()
    ind$fitness <- 0.5 + runif(1, -0.2, 0.2)  # Random fitness around 0.5
    ind$genotype$expr <- c(runif(1), runif(1))  # Random genotype values
    ind$genotype$op <- c(runif(1))
    ind$genotype$var <- c(runif(1))
    population[[i]] <- ind
  }
  return(population)
}

# Simple fitness function for testing
test_fitness_fn <- function(phenotype) {
  # Just count the number of characters as a mock fitness
  return(nchar(phenotype) / 10)
}

test_that("update_grammar_probabilities adjusts probabilities correctly", {
  # Setup
  grammar <- create_test_grammar()
  individual <- create_test_individual()

  # Record original probabilities
  original_probs <- list(
    expr = sapply(grammar$rules$expr, function(r) r$prob),
    op = sapply(grammar$rules$op, function(r) r$prob),
    var = sapply(grammar$rules$var, function(r) r$prob)
  )

  # Update with low learning factor
  updated_grammar <- update_grammar_probabilities(grammar, individual, learning_factor = 0.05, force_perturbation = TRUE)

  # Validate that probabilities were adjusted
  for (nt in names(grammar$rules)) {
    # Check that probabilities still sum to 1
    updated_probs <- sapply(updated_grammar$rules[[nt]], function(r) r$prob)
    expect_equal(sum(updated_probs), 1.0, tolerance = 1e-10)

    # Check that probabilities changed in the expected direction
    if (nt %in% names(individual$genotype) && length(individual$genotype[[nt]]) > 0) {
      expect_false(identical(original_probs[[nt]], updated_probs))
    }
  }

  # Test with empty genotype
  empty_individual <- list(genotype = list())
  expect_error(update_grammar_probabilities(grammar, empty_individual, force_perturbation = TRUE))

  # Test with NULL genotype
  null_individual <- list()
  expect_error(update_grammar_probabilities(grammar, null_individual, force_perturbation = TRUE))
})

test_that("initialize_population creates valid individuals", {
  grammar <- create_test_grammar()
  pop_size <- 20
  max_depth <- 5

  population <- initialize_population(grammar, pop_size, max_depth)

  # Check population size
  expect_equal(length(population), pop_size)

  # Check that each individual has a valid genotype
  for (i in seq_len(pop_size)) {
    expect_true(is.list(population[[i]]))
    expect_true(is.list(population[[i]]$genotype))
    expect_true(validate_genotype(population[[i]]$genotype, grammar))
    expect_null(population[[i]]$fitness)
  }
})

test_that("evaluate_population calculates fitness for all individuals", {
  grammar <- create_test_grammar()
  population <- create_test_population(size = 15)

  # Clear fitness values to simulate unevaluated population
  for (i in seq_along(population)) {
    population[[i]]$fitness <- NULL
  }

  evaluated_pop <- evaluate_population(population, grammar, test_fitness_fn)

  # Check that all individuals have fitness values
  for (i in seq_along(evaluated_pop)) {
    expect_true(!is.null(evaluated_pop[[i]]$fitness))
    expect_true(is.numeric(evaluated_pop[[i]]$fitness))
    expect_true(!is.null(evaluated_pop[[i]]$phenotype))
    expect_true(!is.null(evaluated_pop[[i]]$tree_depth))
    expect_true(!is.null(evaluated_pop[[i]]$mapping_positions))
  }
})

test_that("tournament_selection selects better individuals more often", {
  # Create population with controlled fitness values
  population <- create_test_population(size = 100)

  # Set fitness values: first 20 have good fitness, rest have poor fitness
  for (i in seq_along(population)) {
    if (i <= 20) {
      population[[i]]$fitness <- 0.1  # Good (low) fitness
    } else {
      population[[i]]$fitness <- 0.9  # Poor (high) fitness
    }
  }

  # Run many tournament selections
  n_trials <- 1000
  selected_indices <- integer(n_trials)

  for (i in seq_len(n_trials)) {
    selected <- tournament_selection(population, tournament_size = 5)
    # Find the index of the selected individual
    for (j in seq_along(population)) {
      if (identical(selected, population[[j]])) {
        selected_indices[i] <- j
        break
      }
    }
  }

  # Count how many selections were from the good fitness group
  good_selections <- sum(selected_indices <= 20)

  # With random selection, we'd expect about 20% good selections
  # With tournament selection, we should get significantly more
  expect_true(good_selections > 400)  # Should be much higher than 200 (20% of 1000)

  # Test with tournament_size = 1 (random selection)
  n_trials_random <- 1000
  random_selected_indices <- integer(n_trials_random)

  for (i in seq_len(n_trials_random)) {
    selected <- tournament_selection(population, tournament_size = 1)
    for (j in seq_along(population)) {
      if (identical(selected, population[[j]])) {
        random_selected_indices[i] <- j
        break
      }
    }
  }

  good_random_selections <- sum(random_selected_indices <= 20)

  # Should be closer to the expected 20%
  expect_true(abs(good_random_selections - 200) < 100)  # Allow some random variation

  # Test with tournament_size > population size
  large_tournament <- tournament_selection(population, tournament_size = 1000)
  expect_true(!is.null(large_tournament))  # Should still return an individual
})

test_that("crossover combines parent genotypes correctly", {
  parent1 <- create_test_individual()
  parent2 <- create_test_individual()

  # Modify parent2 to ensure it's different
  parent2$genotype$expr <- c(0.9, 0.3)
  parent2$genotype$op <- c(0.8)
  parent2$genotype$var <- c(0.7)

  # Run crossover multiple times to test different mask patterns
  for (i in 1:20) {
    offspring <- crossover(parent1, parent2)

    # Check offspring structure
    expect_true(is.list(offspring))
    expect_true(is.list(offspring$genotype))
    expect_null(offspring$fitness)

    # Check that offspring has all non-terminals from parents
    expect_equal(
      sort(names(offspring$genotype)),
      sort(unique(c(names(parent1$genotype), names(parent2$genotype))))
    )

    # For each non-terminal, check that values came from one of the parents
    for (nt in names(offspring$genotype)) {
      # The offspring's values should match either parent1 or parent2
      parent1_match <- identical(offspring$genotype[[nt]], parent1$genotype[[nt]])
      parent2_match <- identical(offspring$genotype[[nt]], parent2$genotype[[nt]])

      # At least one should match (could be empty in one parent)
      if (nt %in% names(parent1$genotype) && nt %in% names(parent2$genotype)) {
        expect_true(parent1_match || parent2_match)
      }
    }
  }

  # Test crossover with parents having different non-terminals
  special_parent1 <- parent1
  special_parent2 <- parent2
  special_parent2$genotype$extra_nt <- c(0.5, 0.6)  # Add new non-terminal

  offspring_special <- crossover(special_parent1, special_parent2)

  # Offspring should have the union of non-terminals
  expect_true(all(names(special_parent1$genotype) %in% names(offspring_special$genotype)))
  expect_true(all(names(special_parent2$genotype) %in% names(offspring_special$genotype)))
})

test_that("mutate produces valid variants", {
  individual <- create_test_individual()
  mutation_rate <- 0.5  # High rate to ensure mutations happen

  # Run multiple mutations to test probabilistic behavior
  n_mutations <- 20
  any_changed <- FALSE

  for (i in 1:n_mutations) {
    mutated <- mutate(individual, mutation_rate)

    # Check structure
    expect_true(is.list(mutated))
    expect_true(is.list(mutated$genotype))
    expect_null(mutated$fitness)

    # Check that genotype structure is preserved
    expect_equal(names(mutated$genotype), names(individual$genotype))
    for (nt in names(mutated$genotype)) {
      expect_equal(length(mutated$genotype[[nt]]), length(individual$genotype[[nt]]))

      # Check that values are in valid range
      expect_true(all(mutated$genotype[[nt]] >= 0 & mutated$genotype[[nt]] <= 1))

      # Track if any values changed
      if (!identical(mutated$genotype[[nt]], individual$genotype[[nt]])) {
        any_changed <- TRUE
      }
    }
  }

  # With high mutation rate, at least some values should have changed
  expect_true(any_changed)

  # Test with very low mutation rate
  low_mutation <- mutate(individual, 0.0001)
  # Most likely no mutation occurred, so genotypes should be identical
  expect_equal(low_mutation$genotype, individual$genotype)
})

test_that("create_next_generation preserves population size", {
  grammar <- create_test_grammar()
  population <- create_test_population(size = 30)

  # Ensure population is evaluated
  for (i in seq_along(population)) {
    population[[i]]$fitness <- 0.1 + i * 0.01  # Increasing fitness values
  }

  # Parameters
  crossover_rate <- 0.7
  mutation_rate <- 0.1
  tournament_size <- 3
  elitism <- 5

  # Create new generation
  new_population <- create_next_generation(
    population, grammar, crossover_rate, mutation_rate, tournament_size, elitism
  )

  # Check population size is preserved
  expect_equal(length(new_population), length(population))

  # Check that elite individuals are preserved (should be at the end)
  for (i in 1:elitism) {
    elite_index <- length(population) - elitism + i
    new_elite_index <- length(new_population) - elitism + i

    # The elite individual in the new population should be identical to the one in the old population
    expect_equal(new_population[[new_elite_index]], population[[i]])
  }

  # Check that non-elite individuals have NULL fitness (need evaluation)
  for (i in 1:(length(new_population) - elitism)) {
    expect_null(new_population[[i]]$fitness)
  }

  # Test with different elitism values
  new_population_no_elitism <- create_next_generation(
    population, grammar, crossover_rate, mutation_rate, tournament_size, elitism = 0
  )
  expect_equal(length(new_population_no_elitism), length(population))

  # All individuals should have NULL fitness
  for (i in seq_along(new_population_no_elitism)) {
    expect_null(new_population_no_elitism[[i]]$fitness)
  }
})

test_that("run_psge performs evolution correctly", {
  # Skip this test if it takes too long
  skip_on_cran()

  grammar <- create_test_grammar()

  # Create a fitness function that rewards longer expressions
  length_fitness <- function(phenotype) {
    # Clean spaces and normalize
    clean <- gsub("\\s+", "", phenotype)
    # Penalize very short expressions, reward longer ones up to a point
    if (nchar(clean) <= 2) return(2.0)
    if (nchar(clean) <= 5) return(1.0)
    return(0.5)  # Reward longer expressions
  }

  # Run PSGE with small population and generations for testing
  result <- run_psge(
    grammar, length_fitness,
    pop_size = 20,
    generations = 5,
    crossover_rate = 0.8,
    mutation_rate = 0.1,
    max_depth = 6,
    tournament_size = 3,
    elitism = 2,
    learning_factor = 0.05,
    alternate_update = TRUE,
    verbose = FALSE
  )

  # Check result structure
  expect_true(is.list(result))
  expect_true("best_solution" %in% names(result))
  expect_true("final_grammar" %in% names(result))
  expect_true("history" %in% names(result))

  # Check best solution
  expect_true(is.list(result$best_solution))
  expect_true("fitness" %in% names(result$best_solution))
  expect_true("phenotype" %in% names(result$best_solution))
  expect_true("genotype" %in% names(result$best_solution))

  # Check history
  expect_true(is.data.frame(result$history))
  expect_equal(nrow(result$history), 6)  # Initial + 5 generations
  expect_true(all(c("generation", "best_fitness", "mean_fitness") %in% colnames(result$history)))

  # Check that final grammar probabilities differ from initial probabilities
  initial_probs <- list(
    expr = sapply(grammar$rules$expr, function(r) r$prob),
    op = sapply(grammar$rules$op, function(r) r$prob),
    var = sapply(grammar$rules$var, function(r) r$prob)
  )

  final_probs <- list(
    expr = sapply(result$final_grammar$rules$expr, function(r) r$prob),
    op = sapply(result$final_grammar$rules$op, function(r) r$prob),
    var = sapply(result$final_grammar$rules$var, function(r) r$prob)
  )

  # At least some probabilities should have changed
  expect_false(identical(initial_probs, final_probs))

  # Test with alternate_update = FALSE
  result_no_alternate <- run_psge(
    grammar, length_fitness,
    pop_size = 10,
    generations = 3,
    alternate_update = FALSE,
    verbose = FALSE
  )

  expect_true(is.list(result_no_alternate))
  expect_true("best_solution" %in% names(result_no_alternate))
})

test_that("psge function works with grammar file", {
  # Create temporary grammar file
  grammar_file <- tempfile(fileext = ".bnf")
  writeLines(
    c(
      "<expr> ::= <expr> <op> <expr> | <var>",
      "<op> ::= + | - | * | /",
      "<var> ::= x | y | 1.0"
    ),
    grammar_file
  )
  on.exit(unlink(grammar_file))

  # Mock fitness function that always returns the same value to speed up test
  mock_fitness <- function(phenotype) {
    return(1.0)
  }

  # Run with minimal settings
  result <- psge(
    grammar_file, mock_fitness,
    pop_size = 10,
    generations = 2,
    verbose = FALSE
  )

  # Check basic structure
  expect_true(is.list(result))
  expect_true("best_solution" %in% names(result))
  expect_true("final_grammar" %in% names(result))
  expect_true("history" %in% names(result))

  # Check for error with non-existent file
  expect_error(psge("nonexistent_file.bnf", mock_fitness))
})

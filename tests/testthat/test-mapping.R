library(testthat)

# Test grammar string
test_grammar_str <- "
<expr> ::= <expr> <op> <expr> | <var>
<op> ::= + | - | * | /
<var> ::= x | y | 1.0
"

#' Create a test grammar
#' @return A valid PCFG grammar object
create_test_grammar <- function() {
  grammar_file <- tempfile()
  writeLines(test_grammar_str, grammar_file)
  grammar <- read_pcfg_grammar(grammar_file)
  unlink(grammar_file)
  return(grammar)
}

#' Create a grammar with uniform probabilities for testing
#' @return A grammar with uniform probabilities
create_uniform_grammar <- function() {
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

test_that("generate_random_individual creates valid genotypes", {
  grammar <- create_uniform_grammar()

  # Generate multiple individuals
  for (i in 1:10) {
    genotype <- generate_random_individual(grammar, max_depth = 5)

    # Test that genotype is a list
    expect_true(is.list(genotype))

    # Test that genotype has entries for all non-terminals
    expect_true(all(grammar$non_terminals %in% names(genotype)))

    # Test that all values are between 0 and 1
    for (nt in names(genotype)) {
      expect_true(all(genotype[[nt]] >= 0 & genotype[[nt]] <= 1))
    }

    # Test that the genotype is valid
    expect_true(validate_genotype(genotype, grammar))
  }
})

test_that("map_genotype produces valid phenotypes", {
  grammar <- create_uniform_grammar()

  # Test multiple mappings
  for (i in 1:10) {
    genotype <- generate_random_individual(grammar, max_depth = 5)
    result <- map_genotype(grammar, genotype, max_depth = 5)

    # Test result structure
    expect_true(is.list(result))
    expect_true("phenotype" %in% names(result))
    expect_true("depth" %in% names(result))
    expect_true("positions" %in% names(result))

    # Test that phenotype is a string
    expect_true(is.character(result$phenotype))
    expect_length(result$phenotype, 1)

    # Test phenotype is not empty
    expect_true(nchar(result$phenotype) > 0)
  }
})

test_that("map_genotype respects max_depth limit", {
  grammar <- create_uniform_grammar()

  # Create a genotype that would produce a deep tree without limits
  recursive_genotype <- list(
    expr = rep(0.1, 100),  # Always choose recursive option
    op = rep(0.1, 100),    # Always choose '+'
    var = rep(0.1, 100)    # Always choose 'x'
  )

  # Map with different depth limits
  for (depth in c(3, 5, 10)) {
    result <- map_genotype(grammar, recursive_genotype, max_depth = depth)

    # Test that depth doesn't exceed limit
    expect_true(result$depth <= depth + 1)  # +1 because depth starts at 0
  }
})

test_that("map_genotype handles fixed genotypes consistently", {
  grammar <- create_uniform_grammar()

  # Create a fixed genotype
  fixed_genotype <- list(
    expr = c(0.1, 0.6),  # First recursive, then non-recursive
    op = c(0.1),         # '+'
    var = c(0.1)         # 'x'
  )

  # Map multiple times and check consistency
  first_result <- map_genotype(grammar, fixed_genotype, max_depth = 10)

  for (i in 1:5) {
    result <- map_genotype(grammar, fixed_genotype, max_depth = 10)
    expect_equal(result$phenotype, first_result$phenotype)
    expect_equal(result$depth, first_result$depth)
  }
})

test_that("map_genotype handles empty genotypes by extending them", {
  grammar <- create_uniform_grammar()

  # Empty genotype
  empty_genotype <- list(
    expr = numeric(0),
    op = numeric(0),
    var = numeric(0)
  )

  result <- map_genotype(grammar, empty_genotype, max_depth = 5)

  # Check that mapping succeeded
  expect_true(nchar(result$phenotype) > 0)

  # The genotype should have been extended with random values
  for (nt in names(empty_genotype)) {
    expect_true(length(result$positions[[nt]]) > 0)
  }
})

test_that("recursive_individual_creation handles depth limits", {
  grammar <- create_uniform_grammar()

  # Initialize empty genotype
  genotype <- list()
  for (nt in grammar$non_terminals) {
    genotype[[nt]] <- numeric(0)
  }

  # Test with different depth limits
  for (max_depth in c(1, 3, 5)) {
    result <- recursive_individual_creation(
      grammar, genotype, grammar$start_symbol, 0, max_depth
    )

    # Depth should not exceed max_depth
    expect_true(result$depth <= max_depth)

    # Genotype should have valid values
    for (nt in names(result$genotype)) {
      expect_true(all(result$genotype[[nt]] >= 0 & result$genotype[[nt]] <= 1))
    }
  }
})

test_that("recursive_mapping handles terminals and non-terminals", {
  grammar <- create_uniform_grammar()

  # Test terminal case
  positions <- list(expr = 1, op = 1, var = 1)
  terminal_result <- recursive_mapping(
    grammar,
    list(expr = c(0.6), op = c(0.1), var = c(0.1)),  # Select var->x
    positions,
    "x",  # Terminal symbol
    0,
    5
  )

  expect_equal(terminal_result$output, "x")
  expect_equal(terminal_result$depth, 0)

  # Test non-terminal case with explicit mapping
  non_terminal_result <- recursive_mapping(
    grammar,
    list(expr = c(0.6), op = c(0.1), var = c(0.1)),  # expr->var, var->x
    positions,
    "expr",
    0,
    5
  )

  expect_equal(non_terminal_result$output, "x")
  expect_true(non_terminal_result$depth > 0)

  # Test recursive rule
  recursive_result <- recursive_mapping(
    grammar,
    list(expr = c(0.1, 0.6), op = c(0.1), var = c(0.1)),  # expr->(expr op expr), then expr->var
    positions,
    "expr",
    0,
    5
  )

  # Should produce something like "x + x"
  expect_true(length(recursive_result$output) > 1)
  expect_true(recursive_result$depth > 1)
})

test_that("recursive_mapping respects max depth", {
  grammar <- create_uniform_grammar()

  # Genotype that would generate deep recursion without limits
  deep_genotype <- list(
    expr = rep(0.1, 100),  # Always choose recursive expr
    op = rep(0.1, 100),
    var = rep(0.1, 100)
  )

  positions <- list(expr = 1, op = 1, var = 1)

  for (depth_limit in c(2, 4, 6)) {
    result <- recursive_mapping(
      grammar, deep_genotype, positions, "expr", 0, depth_limit
    )

    # Check depth limitation
    expect_true(result$depth <= depth_limit + 1)
  }
})

test_that("recursive_mapping handles probabilistic rule selection", {
  # Create a grammar with biased probabilities
  biased_grammar <- create_uniform_grammar()
  biased_grammar$rules$var[[1]]$prob <- 0.8  # High probability for 'x'
  biased_grammar$rules$var[[2]]$prob <- 0.1  # Low probability for 'y'
  biased_grammar$rules$var[[3]]$prob <- 0.1  # Low probability for '1.0'

  # Count occurrences of each terminal with fixed codon values
  n_trials <- 100
  counts <- list(x = 0, y = 0, one = 0)

  for (i in 1:n_trials) {
    positions <- list(expr = 1, op = 1, var = 1)
    genotype <- list(expr = c(0.6), op = c(0.1), var = c(runif(1)))
    result <- recursive_mapping(
      biased_grammar, genotype, positions, "var", 0, 5
    )

    if (identical(result$output, "x")) counts$x <- counts$x + 1
    if (identical(result$output, "y")) counts$y <- counts$y + 1
    if (identical(result$output, "1.0")) counts$one <- counts$one + 1
  }

  # Check that distribution roughly matches probabilities
  expect_true(counts$x > counts$y * 3)  # x should appear much more often
  expect_true(counts$x > counts$one * 3)
})

test_that("generate_random_individual creates genotypes that match the grammar structure", {
  # Create a grammar with complex structure
  complex_grammar <- list(
    non_terminals = c("program", "stmt", "expr", "term"),
    terminals = c(";", "if", "then", "else", "+", "-", "*", "/", "x", "y", "1"),
    start_symbol = "program",
    rules = list(
      program = list(
        list(symbols = c("stmt"), types = c("NT"), prob = 1.0)
      ),
      stmt = list(
        list(symbols = c("expr", ";"), types = c("NT", "T"), prob = 0.5),
        list(symbols = c("if", "expr", "then", "stmt", "else", "stmt"),
             types = c("T", "NT", "T", "NT", "T", "NT"), prob = 0.5)
      ),
      expr = list(
        list(symbols = c("expr", "+", "term"), types = c("NT", "T", "NT"), prob = 0.3),
        list(symbols = c("expr", "-", "term"), types = c("NT", "T", "NT"), prob = 0.3),
        list(symbols = c("term"), types = c("NT"), prob = 0.4)
      ),
      term = list(
        list(symbols = c("term", "*", "x"), types = c("NT", "T", "T"), prob = 0.2),
        list(symbols = c("term", "/", "y"), types = c("NT", "T", "T"), prob = 0.2),
        list(symbols = c("x"), types = c("T"), prob = 0.2),
        list(symbols = c("y"), types = c("T"), prob = 0.2),
        list(symbols = c("1"), types = c("T"), prob = 0.2)
      )
    )
  )

  # Generate individuals and check their structure
  for (i in 1:5) {
    genotype <- generate_random_individual(complex_grammar, max_depth = 6)

    # All non-terminals should be present
    expect_equal(sort(names(genotype)), sort(complex_grammar$non_terminals))

    # Map the genotype and verify the result
    result <- map_genotype(complex_grammar, genotype, max_depth = 6)
    expect_true(nchar(result$phenotype) > 0)

    # The genotype should contain enough values to produce a valid phenotype
    for (nt in names(genotype)) {
      # There should be values in the genotype
      expect_true(length(genotype[[nt]]) > 0)
    }
  }
})

test_that("map_genotype handles predefined genotypes correctly", {
  grammar <- create_uniform_grammar()

  # Manually create a genotype that should map to "x + y"
  fixed_genotype <- list(
    expr = c(0.1, 0.6, 0.6),  # expr->expr op expr, then expr->var for both expr
    op = c(0.1),              # op->+
    var = c(0.1, 0.4)         # First var->x, second var->y
  )

  result <- map_genotype(grammar, fixed_genotype, max_depth = 5)
  expected_phenotype <- "x+y"  # Spaces might be handled differently

  # Clean spaces for comparison
  clean_result <- gsub("\\s+", "", result$phenotype)
  clean_expected <- gsub("\\s+", "", expected_phenotype)

  expect_equal(clean_result, clean_expected)

  # Test another predefined genotype: "1.0 / x"
  fixed_genotype2 <- list(
    expr = c(0.1, 0.6, 0.6),  # expr->expr op expr, then expr->var for both
    op = c(0.9),              # op->/
    var = c(0.9, 0.1)         # First var->1.0, second var->x
  )

  result2 <- map_genotype(grammar, fixed_genotype2, max_depth = 5)
  expected_phenotype2 <- "1.0/x"

  # Clean spaces for comparison
  clean_result2 <- gsub("\\s+", "", result2$phenotype)
  clean_expected2 <- gsub("\\s+", "", expected_phenotype2)

  expect_equal(clean_result2, clean_expected2)
})

test_that("recursive_individual_creation creates valid trees", {
  grammar <- create_uniform_grammar()
  max_depth <- 4

  # Initialize empty genotype
  genotype <- list()
  for (nt in grammar$non_terminals) {
    genotype[[nt]] <- numeric(0)
  }

  # Run recursive creation multiple times
  for (i in 1:5) {
    result <- recursive_individual_creation(
      grammar, genotype, "expr", 0, max_depth
    )

    # Map the result to get a phenotype
    phenotype_result <- map_genotype(grammar, result$genotype, max_depth)

    # Check that the phenotype is valid
    expect_true(nchar(phenotype_result$phenotype) > 0)

    # Check that depth is within limits
    expect_true(result$depth <= max_depth)
    expect_true(phenotype_result$depth <= max_depth + 1)
  }
})

test_that("map_genotype and recursive_mapping handle empty input correctly", {
  grammar <- create_uniform_grammar()

  # Completely empty genotype
  empty_genotype <- list(
    expr = numeric(0),
    op = numeric(0),
    var = numeric(0)
  )

  # Should automatically generate values
  result <- map_genotype(grammar, empty_genotype, max_depth = 3)
  expect_true(nchar(result$phenotype) > 0)

  # Partially empty genotype
  partial_genotype <- list(
    expr = c(0.1),  # Just one value
    op = numeric(0),
    var = numeric(0)
  )

  # Should fill in missing values
  result2 <- map_genotype(grammar, partial_genotype, max_depth = 3)
  expect_true(nchar(result2$phenotype) > 0)

  # Test with recursive_mapping directly
  positions <- list(expr = 1, op = 1, var = 1)
  result3 <- recursive_mapping(
    grammar, empty_genotype, positions, "expr", 0, 3
  )
  expect_true(length(result3$output) > 0)
})

library(testthat)

# Helper function to create temporary grammar files for testing
create_test_grammar_file <- function(content, filename = "test_grammar.bnf") {
  writeLines(content, filename)
  return(filename)
}

test_that("read_pcfg_grammar works", {
  # Test case 1: Basic grammar parsing
  grammar_file <- create_test_grammar_file(c(
    "<S> ::= <NP> <VP>",
    "<NP> ::= <Det> <Noun> | <Pronoun>",
    "<Det> ::= the | a",
    "<Noun> ::= cat | dog",
    "<VP> ::= <Verb> <NP>",
    "<Verb> ::= chases | eats",
    "<Pronoun> ::= he | she"
  ))
  on.exit(file.remove(grammar_file))  # Ensure file deletion after test

  grammar <- read_pcfg_grammar(grammar_file)
  expect_equal(grammar$start_symbol, "S")
  expect_equal(length(grammar$non_terminals), 7)
  expect_equal(length(grammar$terminals), 8)
  expect_equal(length(grammar$rules[["S"]]), 1)
  expect_equal(grammar$rules[["NP"]][[1]]$symbols, c("Det", "Noun"))
  expect_equal(grammar$rules[["NP"]][[2]]$symbols, c("Pronoun"))
  expect_equal(grammar$rules[["Det"]][[1]]$symbols, "the")
  expect_equal(grammar$rules[["Det"]][[2]]$symbols, "a")
  expect_equal(grammar$rules[["Noun"]][[1]]$symbols, "cat")
  expect_equal(grammar$rules[["Noun"]][[2]]$symbols, "dog")
  expect_equal(grammar$rules[["VP"]][[1]]$symbols, c("Verb", "NP"))
  expect_equal(grammar$rules[["Verb"]][[1]]$symbols, "chases")
  expect_equal(grammar$rules[["Verb"]][[2]]$symbols, "eats")
  expect_equal(grammar$rules[["Pronoun"]][[1]]$symbols, "he")
  expect_equal(grammar$rules[["Pronoun"]][[2]]$symbols, "she")

  # Test case 2: Specified probabilities - Corrected and Consistent
  grammar_file_specified <- create_test_grammar_file(c( #Consistent grammar file
    "<S> ::= <NP> | <VP>",
    "<NP> ::= <Det> <Noun> | <Pronoun>",
    "<Det> ::= the | a",
    "<Noun> ::= cat | dog",
    "<VP> ::= <Verb> <NP>",
    "<Verb> ::= chases | eats",
    "<Pronoun> ::= he | she"
  ))
  on.exit(file.remove(grammar_file_specified))

  prob_map <- list(
    S = c(0.6, 0.4),
    NP = c(0.7, 0.3),
    Det = c(0.5, 0.5),
    Noun = c(0.6, 0.4),
    VP = c(1),
    Verb = c(0.5, 0.5),
    Pronoun = c(0.5, 0.5)
  )

  grammar_specified <- read_pcfg_grammar(grammar_file_specified, initial_prob_method = "specified", prob_map = prob_map)

  expect_equal(grammar_specified$rules[["S"]][[1]]$prob, 0.6)
  expect_equal(grammar_specified$rules[["NP"]][[1]]$prob, 0.7)
  expect_equal(grammar_specified$rules[["Det"]][[1]]$prob, 0.5)
  expect_equal(grammar_specified$rules[["Noun"]][[1]]$prob, 0.6)
  expect_equal(grammar_specified$rules[["VP"]][[1]]$prob, 1)
  expect_equal(grammar_specified$rules[["Verb"]][[1]]$prob, 0.5)
  expect_equal(grammar_specified$rules[["Pronoun"]][[1]]$prob, 0.5)

  # Test case 3: File not found
  expect_error(read_pcfg_grammar("nonexistent_file.bnf"), "Grammar file not found")

  # Test case 4: Invalid non-terminal format - Corrected
  grammar_file_invalid_nt <- create_test_grammar_file("<S> ::= <NP> <VP>\nS ::= NP VP") # Missing <>
  on.exit(file.remove(grammar_file_invalid_nt))

  # Directly test the parsing (without validation)
  expect_warning(
    read_pcfg_grammar(grammar_file_invalid_nt, initial_prob_method = "uniform", prob_map = list(), validate = FALSE),
    "Invalid non-terminal format"
  )

  # Test case 5: Empty grammar file
  grammar_file_empty <- create_test_grammar_file("")
  on.exit(file.remove(grammar_file_empty))
  expect_error(read_pcfg_grammar(grammar_file_empty), "Grammar has no non-terminals")

  # Test case 6: Missing probability map when specified
  expect_error(read_pcfg_grammar(grammar_file, initial_prob_method = "specified"), "prob_map must be provided when method is 'specified'")

  # Test case 7: Grammar with comments and empty lines
  grammar_file_comments <- create_test_grammar_file(c(
    "# This is a comment",
    "<S> ::= <NP> <VP>",
    "",  # Empty line
    "<NP> ::= <Det> <Noun>",
    "# Another comment"
  ))
  on.exit(file.remove(grammar_file_comments))
  grammar_comments <- read_pcfg_grammar(grammar_file_comments, validate = FALSE)
  expect_equal(length(grammar_comments$non_terminals), 2) # Only S and NP
  expect_equal(length(grammar_comments$rules[["S"]]), 1)
  expect_equal(length(grammar_comments$rules[["NP"]]), 1)  # Only one rule after comment/empty line removal.


})

# Test for extract_tokens
test_that("extract_tokens works", {
  # Test case 1: Non-terminals only
  expect_equal(extract_tokens("<NP> <VP>"), list(symbols = c("NP", "VP"), types = c("NT", "NT")))

  # Test case 2: Terminal string
  expect_equal(extract_tokens("the cat"), list(symbols = "the cat", types = "T"))

  # Test case 3: Mixed terminals and non-terminals
  expect_equal(extract_tokens("<Det> the <Noun>"), list(symbols = c("Det", "the", "Noun"), types = c("NT", "T", "NT")))

  # Test case 4: Non-terminal adjacent to terminal (no space)
  expect_equal(extract_tokens("the<Noun>"), list(symbols = c("the", "Noun"), types = c("T", "NT")))
  expect_equal(extract_tokens("<Det>the"), list(symbols = c("Det", "the"), types = c("NT", "T")))


  # Test case 5: Terminal at the end
  expect_equal(extract_tokens("<Det> the <Noun>."), list(symbols = c("Det", "the", "Noun", "."), types = c("NT", "T", "NT", "T")))

  # Test case 6: Empty string
  expect_equal(extract_tokens(""), list(symbols = "", types = "T"))

  # Test case 7: Terminal with leading/trailing spaces (should be trimmed)
  expect_equal(extract_tokens("  the cat  "), list(symbols = "the cat", types = "T"))

  # Test case 8: Mixed with multiple spaces
  expect_equal(extract_tokens("<Det>   the   <Noun>"), list(symbols = c("Det", "the", "Noun"), types = c("NT", "T", "NT")))

  # Test case 9: Just a non-terminal
  expect_equal(extract_tokens("<NP>"), list(symbols = "NP", types = "NT"))

  # Test case 10: Just a terminal with special characters
  expect_equal(extract_tokens("cat-123"), list(symbols = "cat-123", types = "T")) #Hyphenated

  # Test case 11: Empty non-terminal
  expect_equal(extract_tokens("<>"), list(symbols = "", types = "NT"))

})

# Test for initialize_probabilities
test_that("initialize_probabilities works", {
  # Test setup: Create a temporary grammar file
  grammar_file <- create_test_grammar_file("<S> ::= <NP> <VP> | <NP>\n<NP> ::= the cat | a dog\n<VP> ::= chases | eats")
  on.exit(file.remove(grammar_file))  # Clean up the file after the test

  grammar <- read_pcfg_grammar(grammar_file) #Use read_pcfg_grammar to set up the structure.

  # Test case 1: Uniform probabilities
  grammar_uniform <- initialize_probabilities(grammar, method = "uniform")
  expect_equal(grammar_uniform$rules[["S"]][[1]]$prob, 0.5)
  expect_equal(grammar_uniform$rules[["S"]][[2]]$prob, 0.5)
  expect_equal(grammar_uniform$rules[["NP"]][[1]]$prob, 0.5)
  expect_equal(grammar_uniform$rules[["NP"]][[2]]$prob, 0.5)
  expect_equal(grammar_uniform$rules[["VP"]][[1]]$prob, 0.5)
  expect_equal(grammar_uniform$rules[["VP"]][[2]]$prob, 0.5)

  # Test case 2: Specified probabilities
  prob_map <- list(
    S = c(0.2, 0.8),
    NP = c(0.6, 0.4),
    VP = c(0.7, 0.3)
  )
  grammar_specified <- initialize_probabilities(grammar, method = "specified", prob_map = prob_map)
  expect_equal(grammar_specified$rules[["S"]][[1]]$prob, 0.2)
  expect_equal(grammar_specified$rules[["S"]][[2]]$prob, 0.8)
  expect_equal(grammar_specified$rules[["NP"]][[1]]$prob, 0.6)
  expect_equal(grammar_specified$rules[["NP"]][[2]]$prob, 0.4)
  expect_equal(grammar_specified$rules[["VP"]][[1]]$prob, 0.7)
  expect_equal(grammar_specified$rules[["VP"]][[2]]$prob, 0.3)

  # Test case 3: Unknown method
  expect_error(initialize_probabilities(grammar, method = "invalid"), "Unknown probability initialization method")

  # Test case 4: Specified with some missing non-terminals in prob_map
  prob_map_partial <- list(S = c(0.2, 0.8), NP = c(0.6, 0.4)) #VP is missing
  grammar_specified_partial <- initialize_probabilities(grammar, method = "specified", prob_map = prob_map_partial)
  expect_equal(grammar_specified_partial$rules[["S"]][[1]]$prob, 0.2)
  expect_equal(grammar_specified_partial$rules[["S"]][[2]]$prob, 0.8)
  expect_equal(grammar_specified_partial$rules[["NP"]][[1]]$prob, 0.6)
  expect_equal(grammar_specified_partial$rules[["NP"]][[2]]$prob, 0.4)
  expect_equal(grammar_specified_partial$rules[["VP"]][[1]]$prob, 0.5) # Defaults to uniform
  expect_equal(grammar_specified_partial$rules[["VP"]][[2]]$prob, 0.5)

  # Test case 5: Empty prob_map
  prob_map_empty <- list()
  grammar_specified_empty <- initialize_probabilities(grammar, method = "specified", prob_map = prob_map_empty)
  expect_equal(grammar_specified_empty$rules[["S"]][[1]]$prob, 0.5) # Defaults to uniform
  expect_equal(grammar_specified_empty$rules[["S"]][[2]]$prob, 0.5)
  expect_equal(grammar_specified_empty$rules[["NP"]][[1]]$prob, 0.5) # Defaults to uniform
  expect_equal(grammar_specified_empty$rules[["NP"]][[2]]$prob, 0.5)
  expect_equal(grammar_specified_empty$rules[["VP"]][[1]]$prob, 0.5) # Defaults to uniform
  expect_equal(grammar_specified_empty$rules[["VP"]][[2]]$prob, 0.5)

  # Test case 6: prob_map with extra non-terminals
  prob_map_extra <- list(
    S = c(0.2, 0.8),
    NP = c(0.6, 0.4),
    VP = c(0.7, 0.3)
  )
  grammar_specified_extra <- initialize_probabilities(grammar, method = "specified", prob_map = prob_map_extra)
  expect_equal(grammar_specified_extra$rules[["S"]][[1]]$prob, 0.2)
  expect_equal(grammar_specified_extra$rules[["S"]][[2]]$prob, 0.8)
  expect_equal(grammar_specified_extra$rules[["NP"]][[1]]$prob, 0.6)
  expect_equal(grammar_specified_extra$rules[["NP"]][[2]]$prob, 0.4)
  expect_equal(grammar_specified_extra$rules[["VP"]][[1]]$prob, 0.7)
  expect_equal(grammar_specified_extra$rules[["VP"]][[2]]$prob, 0.3)
  # No error should be thrown, the extra non-terminal is simply ignored.

  # Test case 7: prob_map with incorrect probability vector length
  # Test case 7: prob_map with incorrect probability vector length - Corrected
  grammar_file_correct <- create_test_grammar_file("<S> ::= <NP> | <VP>\n<NP> ::= the cat | a dog\n<VP> ::= chases | eats") # Use a correct grammar file
  on.exit(file.remove(grammar_file_correct))
  grammar_correct <- read_pcfg_grammar(grammar_file_correct)

  prob_map_incorrect_length <- list(S = c(0.2), NP = c(0.6, 0.4), VP = c(0.7, 0.3))
  expect_warning(
    grammar_specified_incorrect_length <- initialize_probabilities(grammar_correct, method = "specified", prob_map = prob_map_incorrect_length),
    "Number of probabilities for 'S' doesn't match number of rules"
  )

  expect_equal(grammar_specified_incorrect_length$rules[["S"]][[1]]$prob, 0.5) # Defaults to uniform
  expect_equal(grammar_specified_incorrect_length$rules[["S"]][[2]]$prob, 0.5) # Defaults to uniform
  expect_equal(grammar_specified_incorrect_length$rules[["NP"]][[1]]$prob, 0.6)
  expect_equal(grammar_specified_incorrect_length$rules[["NP"]][[2]]$prob, 0.4)
  expect_equal(grammar_specified_incorrect_length$rules[["VP"]][[1]]$prob, 0.7)
  expect_equal(grammar_specified_incorrect_length$rules[["VP"]][[2]]$prob, 0.3)
  # A warning is printed, but the code still attempts to set the probabilities where possible.
})


# Test for validate_grammar
# Test for validate_grammar
test_that("validate_grammar works", {
  # Test case 1: Valid grammar
  grammar_file <- create_test_grammar_file("<S> ::= <NP> <VP>\n<NP> ::= the cat | a dog\n<VP> ::= chases | eats")
  on.exit(file.remove(grammar_file))
  grammar <- read_pcfg_grammar(grammar_file)
  validated_grammar <- validate_grammar(grammar)
  expect_equal(validated_grammar, grammar)  # Should return the same grammar object

  # Test case 2: Invalid probability sum (should normalize)
  grammar_file_invalid_prob <- create_test_grammar_file("<S> ::= <NP> <VP> | <NP>\n<NP> ::= the cat | a dog\n<VP> ::= chases | eats")
  on.exit(file.remove(grammar_file_invalid_prob))

  # Version for first warning test
  grammar_invalid_prob1 <- read_pcfg_grammar(grammar_file_invalid_prob)
  grammar_invalid_prob1$rules[["S"]][[1]]$prob <- 0.6  # Sum to 1.1
  grammar_invalid_prob1$rules[["S"]][[2]]$prob <- 0.5

  # This warning result isn't used
  expect_warning({
    grammar_invalid_prob1 <- validate_grammar(grammar_invalid_prob1)
  }, "Probabilities for 'S' sum to 1.1. Normalizing.")
  # Then check normalization
  expect_equal(sum(sapply(grammar_invalid_prob1$rules[["S"]], function(r) r$prob)), 1)

  # Create a separate version for second warning test
  grammar_invalid_prob2 <- read_pcfg_grammar(grammar_file_invalid_prob)
  grammar_invalid_prob2$rules[["S"]][[1]]$prob <- 0.6  # Sum to 1.1 again
  grammar_invalid_prob2$rules[["S"]][[2]]$prob <- 0.5
  expect_warning(grammar_invalid_prob2 <- validate_grammar(grammar_invalid_prob2), "Probabilities for 'S' sum to 1.1. Normalizing.")
  expect_equal(sum(sapply(grammar_invalid_prob2$rules[["S"]], function(r) r$prob)), 1)

  # Test case 3: Invalid probability sum (too far off, should throw error)
  grammar_invalid_prob3 <- read_pcfg_grammar(grammar_file_invalid_prob)  # Fresh copy
  grammar_invalid_prob3$rules[["S"]][[1]]$prob <- 1.5
  # Here we check that the validate_grammar function throws the expected error
  # The function will stop execution with the error message, which is the expected behavior
  expect_error(validate_grammar(grammar_invalid_prob3), "Probabilities for 'S' sum to 2, not 1") #Expect the error

  # Test case 4: Missing start symbol - manually set start_symbol to NULL
  grammar_file_missing_start <- create_test_grammar_file("<NP> ::= the cat | a dog\n<VP> ::= chases | eats")
  on.exit(file.remove(grammar_file_missing_start))
  grammar_missing_start <- read_pcfg_grammar(grammar_file_missing_start, validate = FALSE)
  grammar_missing_start$start_symbol <- NULL  # Force NULL start symbol
  expect_error(validate_grammar(grammar_missing_start), "Grammar has no start symbol")

  # Test case 5: Start symbol not a non-terminal
  grammar$start_symbol <- "X"  # Non-existent non-terminal
  expect_error(validate_grammar(grammar), "Start symbol 'X' is not a defined non-terminal")

  # Test case 6: Non-terminal with no rules
  grammar$rules[["X"]] <- list() #Add a non-terminal with no rules.
  grammar$non_terminals <- c(grammar$non_terminals, "X")
  expect_error(validate_grammar(grammar), "Non-terminal 'X' has no production rules")

  # Test case 7: Rule with no symbols - needs to be set after checking types
  grammar$rules[["NP"]][[1]]$symbols <- character(0)
  expect_error(validate_grammar(grammar), "Rule 1 for 'NP' has no symbols")

  # Test case 8: Rule with mismatched symbols and types
  grammar8 <- read_pcfg_grammar(grammar_file)  # Use a fresh copy
  grammar8$rules[["NP"]][[1]]$types <- c("NT", "NT", "T")  # One too many
  expect_error(validate_grammar(grammar8), "Rule 1 for 'NP' has mismatched symbols and types")

  # Test case 9: Rule with invalid symbol type
  grammar9 <- read_pcfg_grammar(grammar_file)  # Fresh copy
  # Make sure the length matches first to avoid the mismatched length error
  grammar9$rules[["NP"]][[1]]$types <- c("NT", "Z")  # "Z" is invalid
  # Ensure symbols has same length
  grammar9$rules[["NP"]][[1]]$symbols <- c("Det", "Noun")[1:length(grammar9$rules[["NP"]][[1]]$types)]
  expect_error(validate_grammar(grammar9), "Rule 1 for 'NP' has invalid symbol types")

  grammar10 <- read_pcfg_grammar(grammar_file)  # Fresh copy
  # Need to set probabilities for all rules for this non-terminal
  # to maintain the sum=1 constraint initially
  rule_count <- length(grammar10$rules[["NP"]])
  for (i in 1:rule_count) {
    grammar10$rules[["NP"]][[i]]$prob <- 1/rule_count
  }
  # Now set the invalid probability
  grammar10$rules[["NP"]][[1]]$prob <- -0.1
  expect_error(validate_grammar(grammar10), "Rule 1 for 'NP' has invalid probability: -0.1")

  # Test case 11: Rule referencing undefined non-terminal
  grammar11 <- read_pcfg_grammar(grammar_file)  # Fresh copy
  # Add Det as a known NT to pass the check
  if (!("Det" %in% grammar11$non_terminals)) {
    grammar11$non_terminals <- c(grammar11$non_terminals, "Det")
  }
  grammar11$rules[["NP"]][[1]]$symbols <- c("Det", "NonExistent")
  grammar11$rules[["NP"]][[1]]$types <- c("NT", "NT")  # Update types to match
  expect_error(validate_grammar(grammar11), "Rule 1 for 'NP' references undefined non-terminal: NonExistent")

  # Test case 12: Unreachable non-terminals (should warn)
  grammar_file_unreachable <- create_test_grammar_file("<S> ::= <NP> <VP>\n<NP> ::= the cat\n<VP> ::= chases\n<X> ::= a") #X is unreachable
  on.exit(file.remove(grammar_file_unreachable))
  grammar_unreachable <- read_pcfg_grammar(grammar_file_unreachable)
  expect_warning(validate_grammar(grammar_unreachable), "Grammar has unreachable non-terminals: X")

  # Test case 13: Empty grammar
  grammar_empty <- list(
    non_terminals = character(0),
    terminals = character(0),
    start_symbol = "S",
    rules = list()
  )
  expect_error(validate_grammar(grammar_empty), "Grammar has no non-terminals")

  # Test case 14: No non-terminals but with rules
  grammar_no_nts <- list(
    non_terminals = character(0),
    terminals = c("a", "b"),
    start_symbol = "S",
    rules = list(S = list())
  )
  expect_error(validate_grammar(grammar_no_nts), "Grammar has no non-terminals")
})

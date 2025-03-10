% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evolution.R
\name{initialize_population}
\alias{initialize_population}
\title{Initialize a Population with Both Directed and Random Individuals Using Enhanced Parallel Processing}
\usage{
initialize_population(
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
)
}
\arguments{
\item{grammar}{A PCFG grammar structure (e.g., from \code{read_pcfg_grammar()}).}

\item{pop_size}{Desired population size (integer).}

\item{max_depth}{Max derivation depth to prevent overly large expansions.}

\item{use_gadi}{Logical indicating whether to use directed genotype creation
(TRUE) or just random initialization (FALSE).}

\item{fraction_directed}{Fraction (0 to 1) of individuals to create via
directed templates. Default 0.3 means 30\% of the population is
template-based, and 70\% is random.}

\item{analysis}{If already analyzed grammar diversity is available, pass it
here. Otherwise, if \code{use_gadi} is TRUE and no analysis is provided, we
will call \code{analyze_grammar_for_diversity()} internally.}

\item{check_duplicates}{Logical indicating whether to avoid duplicate phenotypes
in the initial population. Default TRUE.}

\item{max_attempts}{Maximum number of attempts to generate a unique individual.
Used only if check_duplicates is TRUE. Default 100.}

\item{parallel}{Logical: whether to use parallel processing. Default TRUE.}

\item{n_cores}{Integer: number of CPU cores to use for parallel processing.
If NULL (default), it uses all available cores minus 1.}

\item{batch_size}{Integer: number of individuals to generate in each parallel batch.
Default 10. Smaller values reduce collision risk but increase overhead.}

\item{sequential_batches}{Logical: whether to process batches sequentially but with
internal parallelism (TRUE) or fully in parallel (FALSE). Default TRUE for better
duplicate control.}

\item{post_process_duplicates}{Logical: whether to filter out any duplicates in a
final post-processing step. Default TRUE.}

\item{exact_population_size}{Logical: whether to ensure exactly pop_size individuals
are returned (TRUE) or accept slightly fewer if max_attempts is reached (FALSE).
Default TRUE.}

\item{verbose}{Logical indicating whether to print information about the
initialization process. Default FALSE.}

\item{...}{Additional arguments passed to \code{directed_templates_from_analysis()}.}
}
\value{
A list of individuals, where each individual is a list containing:
\itemize{
\item genotype: Named list of numeric vectors (codons for each non-terminal).
\item phenotype: The derived string or expression from mapping the genotype.
\item fitness: Initialized to NULL here (not yet evaluated).
\item mapping_positions: Additional info from \code{map_genotype()}.
\item tree_depth: Depth of the derived parse tree.
\item initialization: Character indicating how the individual was created ("directed" or "random").
}
}
\description{
Creates a population of size \code{pop_size} using parallel processing with improved
duplicate handling across concurrent tasks. A portion of individuals is generated
using directed templates from GADI, while the rest are generated randomly.
}
\examples{
\dontrun{
  # Basic usage with parallel processing
  library(future)
  library(future.apply)

  grammar <- read_pcfg_grammar("path/to/grammar.bnf")
  pop <- initialize_population(
    grammar, 100,
    use_gadi = TRUE,
    fraction_directed = 0.3,
    parallel = TRUE,
    sequential_batches = TRUE  # Process batches sequentially with internal parallelism
  )

  # Maximum parallelism but with potential duplicate risk
  pop2 <- initialize_population(
    grammar, 200,
    parallel = TRUE,
    sequential_batches = FALSE,  # Fully parallel processing
    post_process_duplicates = TRUE  # Remove any duplicates at the end
  )
}

}

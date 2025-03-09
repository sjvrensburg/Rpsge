# Rpsge: Probabilistic Structured Grammatical Evolution in R 🧬

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

Rpsge implements Probabilistic Structured Grammatical Evolution (PSGE) in R, providing a full evolutionary framework guided by probabilistic grammars. Beyond basic symbolic regression examples, Rpsge now features a robust suite of tools for:

- **Advanced Grammar Analysis:** Analyze grammar structure with functions like `read_pcfg_grammar()` and `analyze_grammar_for_diversity()` to extract key diversity metrics, rule dependencies, and recursive patterns.
- **Population Initialization:** Create populations using both directed initialization (via Grammar Analysis for Diverse Initialization, GADI) and random individuals. Functions such as `initialize_population()` offer sophisticated duplicate control and parallel generation.
- **Evolutionary Operators:** Apply tournament selection, crossover, and mutation through dedicated functions (e.g., `tournament_selection()`, `crossover()`, and `mutate()`) to generate new generations.
- **Genotype-to-Phenotype Mapping:** Leverage a high-performance C++ backend (via Rcpp) for deterministic mapping of genotypes (codon vectors) to phenotypes, using functions like `map_genotype()` and `recursiveMappingCpp()`.
- **Parallel Evaluation and Progress Reporting:** Fully integrate parallel processing using the `future` and `future.apply` packages along with progress tracking via `progressr` for efficient fitness evaluations.

## Installation

You can install the development version of Rpsge from GitHub:

```r
# install.packages("devtools")
devtools::install_github("sjvrensburg/Rpsge")
```

Ensure that you have the required dependencies installed:
- Rcpp (for C++ integration)
- future, future.apply, and progressr (for parallel processing and progress tracking)

## Basic Usage

Below is an example that demonstrates a full evolutionary cycle, from grammar analysis and population initialization to running the PSGE algorithm with parallel fitness evaluation.

```r
library(Rpsge)
library(progressr)

# Configure progress reporting
handlers(global = TRUE)
handlers("progress")

# Read and analyze a grammar from a BNF file
grammar <- read_pcfg_grammar("path/to/grammar.bnf")
analysis <- analyze_grammar_for_diversity(grammar)

# Initialize a population with a mix of directed (via GADI) and random individuals
pop <- initialize_population(
  grammar = grammar,
  pop_size = 100,
  max_depth = 17,
  use_gadi = TRUE,
  fraction_directed = 0.3,
  parallel = TRUE,
  verbose = TRUE
)

# Define a fitness function for your problem (e.g., symbolic regression)
fitness_fn <- function(phenotype) {
  # Example: Evolve an expression to approximate x^2 + x + 1
  f <- function(x) eval(parse(text = phenotype))
  x <- seq(-1, 1, length.out = 21)
  y_true <- x^2 + x + 1
  y_pred <- sapply(x, f)
  sqrt(mean((y_pred - y_true)^2))
}

# Run the PSGE algorithm with parallel evaluation and progress reporting
result <- run_psge(
  grammar = grammar,
  fitness_fn = fitness_fn,
  pop_size = 100,
  generations = 50,
  crossover_rate = 0.9,
  mutation_rate = 0.1,
  max_depth = 17,
  tournament_size = 3,
  learning_factor = 0.01,
  force_perturbation = FALSE,
  verbose = TRUE,
  parallel = TRUE,
  n_workers = 2
)

# View the best solution and final grammar
print(result$best_solution$phenotype)
```

## Features

- **Probabilistic Grammar-Guided Evolution:** Utilize a probabilistic context-free grammar (PCFG) to steer the evolutionary search.
- **Advanced Grammar Analysis:** Identify recursive non-terminals, compute rule diversity indices, and determine dependencies using dedicated analysis routines.
- **Directed and Random Initialization:** Generate individuals using both directed templates (GADI) and random sampling, with robust duplicate control.
- **Evolutionary Operators:** Use tournament selection, crossover, and mutation functions to evolve your population.
- **High-Performance C++ Backend:** Benefit from efficient genotype-to-phenotype mapping and rule evaluation powered by Rcpp and Armadillo.
- **Parallel Processing:** Leverage parallel evaluation and batch processing via the `future` and `future.apply` packages, with real-time progress reporting using `progressr`.

## Citation

If you use Rpsge in your research, please cite:

```bibtex
@inproceedings{megane2022psge,
  title={Probabilistic Structured Grammatical Evolution},
  author={M{\'e}gane, Jessica and Louren{\c{c}}o, Nuno and Machado, Penousal},
  booktitle={2022 IEEE Congress on Evolutionary Computation (CEC)},
  pages={1--9},
  year={2022},
  doi={10.1109/CEC55065.2022.9870397}
}
```

## Contributing

Contributions are welcome! Please submit Pull Requests or open issues on GitHub. For developers, please refer to the in-code documentation and tests for guidance on extending the package.

## License

GPL (>= 3)

---

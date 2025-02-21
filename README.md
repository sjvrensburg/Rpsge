# Rpsge: Probabilistic Structured Grammatical Evolution in R 🧬

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview 🔍

Rpsge implements Probabilistic Structured Grammatical Evolution (PSGE) in R. PSGE combines the structured representation from Structured Grammatical Evolution with probabilistic grammar learning, creating an adaptive approach to grammar-guided genetic programming.

## Installation 📦

You can install the development version of Rpsge from GitHub:

```r
# install.packages("devtools")
devtools::install_github("sjvrensburg/Rpsge")
```

## Basic Usage 💻

Here's a simple example of solving a symbolic regression problem:

```r
library(Rpsge)
library(progressr)

# Configure progress reporting
handlers(global = TRUE)
handlers("progress")

# Define fitness function
fitness_fn <- function(phenotype) {
  # Example: Try to evolve x^2 + x + 1
  f <- function(x) eval(parse(text = phenotype))
  x <- seq(-1, 1, by = 0.1)
  y_true <- x^2 + x + 1
  y_pred <- sapply(x, f)
  sqrt(mean((y_pred - y_true)^2))
}

# Read grammar file
grammar <- read_pcfg_grammar("path/to/grammar.bnf")

# Run PSGE with progress reporting
with_progress({
  result <- psge(
    grammar_file = "path/to/grammar.bnf",
    fitness_fn = fitness_fn,
    pop_size = 100,
    generations = 50
  )
})

# Print best solution
print(result$best_solution$phenotype)
```

## Features ⭐

- Probabilistic grammar-guided evolution
- High-performance C++ core through Rcpp
- Parallel fitness evaluation support
- Comprehensive grammar validation
- Progress tracking and visualization

## Citation 📚

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

## Contributing 🤝

Contributions are welcome! Please feel free to submit a Pull Request.

## License 📄

GPL (>= 3)

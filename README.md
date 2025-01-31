# Rpsge: Probabilistic Structured Grammatical Evolution in R

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

An R package implementing **Probabilistic Structured Grammatical Evolution (PSGE)**, combining grammatical evolution with probabilistic context-free grammars (PCFG) to evolve solutions for optimization problems.


> :warning: This package still requires extensive testing. Consider it **pre-alpha software**. Use at own risk!

## Features

- 🧬 **Probabilistic Grammar-Guided Evolution**: Uses PCFGs to bias search toward promising regions
- 🔄 **Adaptive Learning**: Dynamically updates grammar probabilities based on successful individuals
- 🏗️ **Structured Genotype-Phenotype Mapping**: Maintains tree structure during evolution
- 🛠️ **Customizable Operators**: Tournament selection, Gaussian mutation, masked crossover
- 📜 **BNF Grammar Support**: Read standard Backus-Naur Form (BNF) grammars
- 🚀 **C++ Accelerated Core**: Critical probability updates optimized with RcppArmadillo

## Important Requirements & Caveats

- Requires a C++ compiler (Rtools on Windows)
- Memory usage scales with population size and grammar complexity
- Grammar files must use UTF-8 encoding
- Solutions evolve stochastically - multiple runs recommended for best results
- Large grammars may require significant computation time

## Installation

```r
# Install dependencies first
install.packages(c("Rcpp", "RcppArmadillo", "R6"))

# Then install Rpsge
if (!require("devtools")) install.packages("devtools")
devtools::install_github("sjvrensburg/Rpsge")
```

## Usage

### 1. Define a Grammar
Create a BNF grammar file (`example.bnf`):
```bnf
<expr> ::= <var> | ( <expr> <op> <expr> )
<op>   ::= + | - | *
<var>  ::= x | 1.0 | 2.0
```

### 2. Create PSGE Instance
```r
library(Rpsge)

# Initialize grammar
grammar <- create_grammar("example.bnf")

# Evaluation function (minimize MSE from f(x) = x + 1)
eval_func <- function(phenotype) {
  tryCatch({
    x <- seq(-1, 1, by = 0.1)
    y_true <- x + 1
    y_pred <- eval(parse(text = phenotype))
    mean((y_true - y_pred)^2)
  }, error = function(e) {
    return(Inf)  # Return worst fitness on error
  })
}

# Create PSGE optimizer
psge <- create_psge(
  grammar = grammar,
  eval_func = eval_func,
  population_size = 100,
  generations = 50,
  learning_factor = 0.1
)
```

### 3. Run Evolution
```r
result <- psge$evolve(verbose = TRUE)
cat("Best solution:", result$phenotype, "\nFitness:", result$fitness)
```

## Key Components

### `PSGE` R6 Class
- `$evolve()`: Run evolutionary process
- `$get_best()`: Retrieve best solution
- `$get_population()`: Access current population

### `Grammar` R6 Class
- `$read_grammar()`: Load BNF grammar
- `$generate_uniform_pcfg()`: Initialize equal rule probabilities
- `$mapping()`: Convert genotype to phenotype
- `$find_shortest_path()`: Calculate non-recursive expansions

## Advanced Configuration

### Custom Evolutionary Parameters
```r
psge <- create_psge(
  grammar = grammar,
  eval_func = eval_func,
  population_size = 200,
  generations = 100,
  elite_size = 20,
  tournament_size = 5,
  crossover_prob = 0.85,
  mutation_prob = 0.15,
  learning_factor = 0.2
)
```

### Grammar Probabilities
Access and modify PCFG directly:
```r
# Get probabilities for <expr> non-terminal
expr_probs <- grammar$pcfg[1, ]

# Force uniform redistribution
grammar$generate_uniform_pcfg()
```

## Documentation

Full documentation available via:
```r
?PSGE
?Grammar
?create_psge
?create_grammar
```

## Dependencies

- R (≥ 4.2.2)
- Rcpp
- RcppArmadillo
- R6

## Contributing

Contributions are welcome! Please submit issues and pull requests through GitHub:
1. Fork repository
2. Create feature branch (`git checkout -b feature/your-feature`)
3. Commit changes (`git commit -am 'Add awesome feature'`)
4. Push to branch (`git push origin feature/your-feature`)
5. Open Pull Request

## License
GNU GPLv3 © Stéfan Janse van Rensburg

## Acknowledgments
This package is based on the Probabilistic Structured Grammatical Evolution (PSGE) algorithm introduced in the following work:

```bibtex
@inproceedings{megane2022cec,
    author={Mégane, Jessica and Lourenço, Nuno and Machado, Penousal},
    booktitle={2022 IEEE Congress on Evolutionary Computation (CEC)}, 
    title={Probabilistic Structured Grammatical Evolution}, 
    year={2022},
    pages={1-9},
    doi={10.1109/CEC55065.2022.9870397}}
```
The implementation in this package is heavily based on Jessica Mégane's original Python code, available at:
[https://github.com/jessicamegane/psge](https://github.com/jessicamegane/psge)
For production use, we highly recommend that you rather use their Python package.

We thank the authors for their work and for making their code publicly available.

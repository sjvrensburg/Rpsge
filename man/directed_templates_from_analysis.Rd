% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/directed_genotypes.R
\name{directed_templates_from_analysis}
\alias{directed_templates_from_analysis}
\title{Create Directed Genotype Templates from Grammar Analysis}
\usage{
directed_templates_from_analysis(
  analysis,
  max_diverse_nt = 5,
  max_rules_per_nt = 3,
  include_start_symbol = TRUE,
  prioritize_uncommon = TRUE,
  create_combined_templates = TRUE,
  dedup_templates = TRUE
)
}
\arguments{
\item{analysis}{The list returned by analyze_grammar_for_diversity().}

\item{max_diverse_nt}{Integer: how many non-terminals to target (in
order of descending diversity). Defaults to 5.}

\item{max_rules_per_nt}{Integer: maximum number of distinct rules
to create templates for in each targeted non-terminal. Defaults to 3.}

\item{include_start_symbol}{Logical: whether to always create templates
that select each rule of the start symbol. Defaults to TRUE.}

\item{prioritize_uncommon}{Logical: whether to prioritize less common rules
(lower probability rules) when selecting templates. Defaults to TRUE.}

\item{create_combined_templates}{Logical: whether to create templates that
combine choices for multiple non-terminals. Defaults to TRUE.}

\item{dedup_templates}{Logical: whether to remove duplicate templates.
Defaults to TRUE.}
}
\value{
A list of partial templates. Each element is itself a named list,
e.g. \code{list(expr = 0.3, product_expr = 0.6, ...)} indicating
codons for some non-terminals. Others remain \code{NULL} (or omitted),
allowing for random filling later.
}
\description{
Produces a set of partial genotype "templates" that bias
the GE mapping process toward high-diversity grammar paths,
with particular emphasis on less common rules.
}

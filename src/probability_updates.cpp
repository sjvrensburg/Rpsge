#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat update_pcfg_probabilities(arma::mat& pcfg,
                               const arma::mat& mask,
                               const std::vector<std::vector<int>>& grammar_counter,
                               double learning_factor) {
  int rows = pcfg.n_rows;
  int cols = pcfg.n_cols;

  // For each non-terminal
  for(int i = 0; i < rows; i++) {
    // Skip if only one production rule
    if(arma::accu(mask.row(i)) <= 1) continue;

    // Get total expansions for this non-terminal
    double total = 0;
    for(int count : grammar_counter[i]) {
      total += count;
    }

    if(total == 0) continue;

    // Update probabilities for each production rule
    for(int j = 0; j < cols; j++) {
      if(!mask(i,j)) continue;

      double counter = grammar_counter[i][j];
      double old_prob = pcfg(i,j);

      if(counter > 0) {
        // Rule was used - increase probability
        pcfg(i,j) = std::min(old_prob + learning_factor * counter / total, 1.0);
      } else {
        // Rule was not used - decrease probability
        pcfg(i,j) = std::max(old_prob - learning_factor * old_prob, 0.0);
      }
    }

    // Normalize row probabilities to sum to 1
    double row_sum = arma::accu(arma::clamp(pcfg.row(i), 0, arma::datum::inf));
    if(row_sum > 0) {
      pcfg.row(i) = arma::clamp(pcfg.row(i), 0, arma::datum::inf) / row_sum;
    }
  }
  return pcfg;
}

// [[Rcpp::export]]
std::vector<std::vector<int>> get_grammar_counter(const List& individual,
                                                  const CharacterVector& non_terminals,
                                                  const List& grammar) {
  int n_nts = non_terminals.length();
  std::vector<std::vector<int>> counter(n_nts);

  // Initialize counters for each non-terminal
  for(int i = 0; i < n_nts; i++) {
    String nt = non_terminals[i];
    List productions = grammar[nt];
    counter[i] = std::vector<int>(productions.length(), 0);
  }

  // Count expansions from genotype
  List genotype = individual["genotype"];
  for(int i = 0; i < n_nts; i++) {
    String nt = non_terminals[i];
    List expansions = genotype[nt];

    for(int j = 0; j < expansions.length(); j++) {
      NumericVector expansion = expansions[j];
      counter[i][expansion[0]]++;
    }
  }

  return counter;
}

// [[Rcpp::export]]
List independent_update(List& grammar_obj,
                        const List& best,
                        double learning_factor) {

  CharacterVector non_terminals = grammar_obj.attr("non_terminals");
  List grammar = grammar_obj.attr("grammar");
  arma::mat pcfg = grammar_obj.attr("pcfg");
  arma::mat mask = grammar_obj.attr("pcfg_mask");

  // Get expansion counts from best individual
  std::vector<std::vector<int>> counter = get_grammar_counter(best, non_terminals, grammar);

  // Update PCFG probabilities
  update_pcfg_probabilities(pcfg, mask, counter, learning_factor);

  // Update grammar object attributes
  grammar_obj.attr("pcfg") = pcfg;

  return grammar_obj;
}

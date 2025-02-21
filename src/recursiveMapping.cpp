#include <Rcpp.h>
#include <random>
#include <string>
#include <sstream>
#include <vector>

using namespace Rcpp;

// Helper: compute deterministic seed from position information
std::uint32_t computePosHash(const List& positions, const std::string& symbol, int current_depth) {
  std::ostringstream oss;
  CharacterVector posNames = positions.attr("names");

  // Build hash string exactly as R version does
  for (int i = 0; i < positions.length(); i++) {
    std::string name = Rcpp::as<std::string>(posNames[i]);
    double pos = Rcpp::as<double>(positions[i]);
    oss << name << ":" << pos << ";";
  }
  oss << "NT:" << symbol << ";depth:" << current_depth;

  // Sum UTF-8 codes as in R
  std::string hashStr = oss.str();
  std::uint32_t sum = 0;
  for (unsigned char ch : hashStr) {
    sum += ch;
  }
  return sum;
}

// Helper: generate random number with optional seeding
double generateRandom(double seed, bool useSeed) {
  if (useSeed) {
    std::mt19937 gen(static_cast<unsigned int>(seed));
    std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
  } else {
    thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    return dis(gen);
  }
}

// Helper: check if a rule is recursive with respect to the current non-terminal
bool isRuleRecursive(const List& rule, const std::string& non_terminal) {
  CharacterVector symbols = rule["symbols"];
  CharacterVector types = rule["types"];

  for (int i = 0; i < symbols.length(); i++) {
    if (as<std::string>(types[i]) == "NT" &&
        as<std::string>(symbols[i]) == non_terminal) {
      return true;
    }
  }
  return false;
}

// Helper: select rule based on codon and probabilities
List selectRule(const List& rules, double codon, bool debug_output = false) {
  int n = rules.length();
  std::vector<double> probs(n);
  double total = 0.0;

  // Extract probabilities
  for (int i = 0; i < n; i++) {
    List rule = rules[i];
    double p = as<double>(rule["prob"]);
    probs[i] = p;
    total += p;
  }

  // Handle zero probabilities case
  if (total == 0.0) {
    std::fill(probs.begin(), probs.end(), 1.0 / n);
    total = 1.0;
  }

  if (debug_output) {
    Rcpp::Rcout << "Codon: " << codon << ", Rules: " << n << ", Total prob: " << total << "\n";
    for (int i = 0; i < n; i++) {
      List rule = rules[i];
      CharacterVector syms = rule["symbols"];
      Rcpp::Rcout << "Rule " << i << ": ";
      for (int j = 0; j < syms.length(); j++) {
        Rcpp::Rcout << as<std::string>(syms[j]) << " ";
      }
      Rcpp::Rcout << "prob=" << probs[i] << "\n";
    }
  }

  // Select rule
  double cum = 0.0;
  int selected = n - 1;  // Default to last rule

  if (debug_output) {
    Rcpp::Rcout << "Cumulative thresholds: ";
  }

  for (int i = 0; i < n; i++) {
    cum += probs[i] / total;
    if (debug_output) {
      Rcpp::Rcout << cum << " ";
    }
    if (codon <= cum) {
      selected = i;
      break;
    }
  }

  if (debug_output) {
    Rcpp::Rcout << "\nSelected rule " << selected << "\n";
    List sel = rules[selected];
    CharacterVector syms = sel["symbols"];
    Rcpp::Rcout << "Selected symbols: ";
    for (int i = 0; i < syms.length(); i++) {
      Rcpp::Rcout << as<std::string>(syms[i]) << " ";
    }
    Rcpp::Rcout << "\n";
  }

  return rules[selected];
}

// Main recursive mapping function
// [[Rcpp::export]]
List recursiveMappingCpp(List grammar, List genotype, List positions,
                         std::string symbol, int current_depth, int max_depth,
                         CharacterVector output, bool use_position_seeds = false,
                         bool debug_output = false) {

  if (debug_output) {
    Rcpp::Rcout << "Mapping symbol: " << symbol
                << " at depth: " << current_depth << "\n";
  }

  // Guard against exceeding max depth
  if (current_depth > max_depth) {
    output.push_back("max_depth_exceeded");
    return List::create(
      _["output"] = output,
      _["depth"] = current_depth,
      _["positions"] = positions,
      _["genotype"] = genotype
    );
  }

  // Check if symbol is a non-terminal
  CharacterVector nonTerms = grammar["non_terminals"];
  bool isNonTerminal = false;
  for (int i = 0; i < nonTerms.length(); i++) {
    if (as<std::string>(nonTerms[i]) == symbol) {
      isNonTerminal = true;
      break;
    }
  }

  // If not a non-terminal, add to output and return
  if (!isNonTerminal) {
    output.push_back(symbol);
    return List::create(
      _["output"] = output,
      _["depth"] = current_depth,
      _["positions"] = positions,
      _["genotype"] = genotype
    );
  }

  // Get current position and check genotype
  double pos = as<double>(positions[symbol]);
  NumericVector genVec = genotype[symbol];
  int pos_idx = static_cast<int>(pos) - 1;  // Convert to 0-based index

  if (debug_output) {
    Rcpp::Rcout << "Position for " << symbol << ": " << pos
                << " (index " << pos_idx << ")\n";
    Rcpp::Rcout << "Available codons: ";
    for (int i = 0; i < genVec.length(); i++) {
      Rcpp::Rcout << genVec[i] << " ";
    }
    Rcpp::Rcout << "\n";
  }

  // Generate new codon if needed
  if (pos_idx >= genVec.length()) {
    double new_value;
    if (use_position_seeds) {
      std::uint32_t seed = computePosHash(positions, symbol, current_depth);
      new_value = generateRandom(seed, true);
    } else {
      new_value = generateRandom(0, false);
    }
    genVec.push_back(new_value);
    genotype[symbol] = genVec;
  }

  // Get codon and update position
  double codon = genVec[pos_idx];
  positions[symbol] = pos + 1;

  if (debug_output) {
    Rcpp::Rcout << "Using codon at position " << pos_idx
                << ": " << codon << "\n";
  }
  // Handle NA/NaN codon
  if (ISNA(codon) || ISNAN(codon)) {
    if (use_position_seeds) {
      // Get position index (1-based to match R)
      CharacterVector posNames = positions.names();
      int pos_index = 0;
      for (int i = 0; i < posNames.length(); i++) {
        if (as<std::string>(posNames[i]) == symbol) {
          pos_index = i + 1;
          break;
        }
      }
      double seed = pos * 1000 + pos_index * 100;
      codon = generateRandom(seed, true);
    } else {
      codon = generateRandom(0, false);
    }
    genVec[static_cast<int>(pos) - 1] = codon;
    genotype[symbol] = genVec;
  }

  // Clamp codon to [0,1]
  codon = std::max(0.0, std::min(1.0, codon));

  if (debug_output) {
    Rcpp::Rcout << "\nSelecting rule for symbol: " << symbol
                << " at depth: " << current_depth
                << " with codon: " << codon << "\n";
  }

  // Get rules and select one
  List allRules = grammar["rules"];
  if (debug_output) {
    Rcpp::Rcout << "Grammar rules structure:\n";
    for (int i = 0; i < allRules.length(); i++) {
      Rcpp::CharacterVector allRulesNames = allRules.names();
      std::string name = "";
      if (allRulesNames.size() != 0) name = Rcpp::as<std::string>(allRulesNames[i]);
      Rcpp::Rcout << "  " << name << ":\n";
      List rules = allRules[i];
      for (int j = 0; j < rules.length(); j++) {
        List rule = rules[j];
        CharacterVector syms = rule["symbols"];
        double p = as<double>(rule["prob"]);
        Rcpp::Rcout << "    Rule " << j << " (p=" << p << "): ";
        for (int k = 0; k < syms.length(); k++) {
          Rcpp::Rcout << as<std::string>(syms[k]) << " ";
        }
        Rcpp::Rcout << "\n";
      }
    }
  }
  List rulesForSymbol = allRules[symbol];
  List selected_rule;

  // At max depth, use non-recursive rules if available
  if (current_depth >= max_depth) {
    std::vector<List> nonRecRules;
    for (int i = 0; i < rulesForSymbol.length(); i++) {
      List rule = rulesForSymbol[i];
      if (!isRuleRecursive(rule, symbol)) {
        nonRecRules.push_back(rule);
      }
    }

    // If no non-recursive rules, use all rules
    if (nonRecRules.empty()) {
      selected_rule = selectRule(rulesForSymbol, codon, debug_output);
    } else {
      List nonRecList(nonRecRules.begin(), nonRecRules.end());
      selected_rule = selectRule(nonRecList, codon, debug_output);
    }
  } else {
    selected_rule = selectRule(rulesForSymbol, codon, debug_output);
  }

  // Check if rule is recursive
  bool is_recursive_rule = isRuleRecursive(selected_rule, symbol);

  // Process selected rule's symbols
  CharacterVector ruleSymbols = selected_rule["symbols"];
  CharacterVector ruleTypes = selected_rule["types"];
  int max_depth_reached = current_depth;

  for (int i = 0; i < ruleSymbols.length(); i++) {
    std::string sym = as<std::string>(ruleSymbols[i]);
    std::string typ = as<std::string>(ruleTypes[i]);

    if (typ == "NT") {
      if (current_depth + 1 <= max_depth) {
        List res = recursiveMappingCpp(
          grammar, genotype, positions, sym,
          current_depth + 1, max_depth, output,
          use_position_seeds, debug_output
        );
        output = res["output"];
        positions = res["positions"];
        genotype = res["genotype"];
        max_depth_reached = std::max(
          max_depth_reached,
          as<int>(res["depth"])
        );
      } else {
        output.push_back("at_max_depth");
        max_depth_reached = max_depth;
      }
    } else {
      output.push_back(sym);
    }
  }

  // Update depth for recursive rules
  if (is_recursive_rule) {
    max_depth_reached = std::max(max_depth_reached, current_depth + 1);
  }

  return List::create(
    _["output"] = output,
    _["depth"] = max_depth_reached,
    _["positions"] = positions,
    _["genotype"] = genotype
  );
}

// Debugging function
// [[Rcpp::export]]
List traceRecursiveMappingWithCodon(List grammar,
                                    NumericVector expr_codons,
                                    NumericVector op_codons,
                                    NumericVector var_codons) {
  List genotype = List::create(
    _["expr"] = expr_codons,
    _["op"] = op_codons,
    _["var"] = var_codons
  );

  List positions = List::create(
    _["expr"] = 1.0,
    _["op"] = 1.0,
    _["var"] = 1.0
  );

  CharacterVector output;
  return recursiveMappingCpp(
    grammar, genotype, positions, "expr", 0, 5,
    output, false, true
  );
}

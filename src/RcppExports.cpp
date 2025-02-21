// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// recursiveMappingCpp
List recursiveMappingCpp(List grammar, List genotype, List positions, std::string symbol, int current_depth, int max_depth, CharacterVector output, bool use_position_seeds, bool debug_output);
RcppExport SEXP _Rpsge_recursiveMappingCpp(SEXP grammarSEXP, SEXP genotypeSEXP, SEXP positionsSEXP, SEXP symbolSEXP, SEXP current_depthSEXP, SEXP max_depthSEXP, SEXP outputSEXP, SEXP use_position_seedsSEXP, SEXP debug_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type grammar(grammarSEXP);
    Rcpp::traits::input_parameter< List >::type genotype(genotypeSEXP);
    Rcpp::traits::input_parameter< List >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< std::string >::type symbol(symbolSEXP);
    Rcpp::traits::input_parameter< int >::type current_depth(current_depthSEXP);
    Rcpp::traits::input_parameter< int >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type output(outputSEXP);
    Rcpp::traits::input_parameter< bool >::type use_position_seeds(use_position_seedsSEXP);
    Rcpp::traits::input_parameter< bool >::type debug_output(debug_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(recursiveMappingCpp(grammar, genotype, positions, symbol, current_depth, max_depth, output, use_position_seeds, debug_output));
    return rcpp_result_gen;
END_RCPP
}
// traceRecursiveMappingWithCodon
List traceRecursiveMappingWithCodon(List grammar, NumericVector expr_codons, NumericVector op_codons, NumericVector var_codons);
RcppExport SEXP _Rpsge_traceRecursiveMappingWithCodon(SEXP grammarSEXP, SEXP expr_codonsSEXP, SEXP op_codonsSEXP, SEXP var_codonsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type grammar(grammarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type expr_codons(expr_codonsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type op_codons(op_codonsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type var_codons(var_codonsSEXP);
    rcpp_result_gen = Rcpp::wrap(traceRecursiveMappingWithCodon(grammar, expr_codons, op_codons, var_codons));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rpsge_recursiveMappingCpp", (DL_FUNC) &_Rpsge_recursiveMappingCpp, 9},
    {"_Rpsge_traceRecursiveMappingWithCodon", (DL_FUNC) &_Rpsge_traceRecursiveMappingWithCodon, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rpsge(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

#pragma once

#include "types/expression_data.h"
#include "types/gene_set.h"
#include <vector>
#include <span>

using namespace std;

namespace gsea {

[[nodiscard]] vector<size_t> generate_random_gene_rank(
    const ExpressionData& expression,
    size_t disease_size);

[[nodiscard]] vector<vector<double>> compute_null_distribution(
    const ExpressionData& expression,
    span<const GeneSet> gene_sets,
    size_t disease_size,
    size_t sample_size);

[[nodiscard]] vector<size_t> find_significant_sets(
    span<const double> actual_scores,
    span<const vector<double>> null_distribution,
    double p_value,
    size_t num_sets);

} // namespace gsea
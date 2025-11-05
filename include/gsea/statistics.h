#pragma once

#include "types/expression_data.h"
#include "types/gene_set.h"
#include <vector>
#include <span>

namespace gsea {

[[nodiscard]] std::vector<size_t> generate_random_gene_rank(
    const ExpressionData& expression,
    size_t disease_size);

[[nodiscard]] std::vector<std::vector<double>> compute_null_distribution(
    const ExpressionData& expression,
    std::span<const GeneSet> gene_sets,
    size_t disease_size,
    size_t sample_size);

[[nodiscard]] std::vector<size_t> find_significant_sets(
    std::span<const double> actual_scores,
    std::span<const std::vector<double>> null_distribution,
    double p_value,
    size_t num_sets);

} // namespace gsea
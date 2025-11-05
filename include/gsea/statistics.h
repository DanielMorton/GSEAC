#pragma once

#include "types/expression_data.h"
#include "types/gene_set.h"
#include <vector>

namespace gsea {

    std::vector<size_t> generate_random_gene_rank(const ExpressionData& expression,
                                                   size_t disease_size);

    std::vector<std::vector<double>> compute_null_distribution(
        const ExpressionData& expression,
        const std::vector<GeneSet>& gene_sets,
        size_t disease_size,
        size_t sample_size);

    std::vector<size_t> find_significant_sets(
        const std::vector<double>& actual_scores,
        const std::vector<std::vector<double>>& null_distribution,
        double p_value,
        size_t num_sets);

} // namespace gsea
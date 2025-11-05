#include "gsea/ranking.h"
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace gsea {

static double calculate_mean(const ExpressionData& expression,
                             size_t gene_idx,
                             std::span<const size_t> sample_indices) {
    double sum = std::accumulate(sample_indices.begin(), sample_indices.end(), 0.0,
        [&](double acc, size_t col) {
            return acc + expression.values()(gene_idx, col);
        });
    return sum / sample_indices.size();
}

std::vector<size_t> compute_gene_rank(const ExpressionData& expression,
                                       std::span<const size_t> disease_indices,
                                       std::span<const size_t> healthy_indices) {
    if (disease_indices.empty() || healthy_indices.empty()) {
        throw std::invalid_argument("Cannot compute gene rank with empty sample groups");
    }

    size_t num_genes = expression.num_genes();

    // Calculate differential expression for each gene
    std::vector<std::pair<size_t, double>> gene_diffs;
    gene_diffs.reserve(num_genes);

    for (size_t gene_idx = 0; gene_idx < num_genes; ++gene_idx) {
        double disease_mean = calculate_mean(expression, gene_idx, disease_indices);
        double healthy_mean = calculate_mean(expression, gene_idx, healthy_indices);
        gene_diffs.emplace_back(gene_idx, disease_mean - healthy_mean);
    }

    // Sort by difference in descending order
    std::ranges::sort(gene_diffs, std::ranges::greater{}, &std::pair<size_t, double>::second);

    // Extract indices
    std::vector<size_t> ranked_indices;
    ranked_indices.reserve(num_genes);
    for (const auto& [idx, _] : gene_diffs) {
        ranked_indices.push_back(idx);
    }

    return ranked_indices;
}

} // namespace gsea
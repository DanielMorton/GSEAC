#include "gsea/ranking.h"
#include <algorithm>
#include <stdexcept>
#include <ranges>
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
    auto gene_diffs = std::views::iota(size_t{0}, num_genes)
        | std::views::transform([&](size_t gene_idx) {
            double disease_mean = calculate_mean(expression, gene_idx, disease_indices);
            double healthy_mean = calculate_mean(expression, gene_idx, healthy_indices);
            return std::pair{gene_idx, disease_mean - healthy_mean};
        })
        | std::ranges::to<std::vector>();

    // Sort by difference in descending order
    std::ranges::sort(gene_diffs, std::ranges::greater{}, &std::pair<size_t, double>::second);

    // Extract indices
    return gene_diffs
        | std::views::elements<0>
        | std::ranges::to<std::vector>();
}

} // namespace gsea
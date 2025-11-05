#include "gsea/ranking.h"
#include <algorithm>
#include <stdexcept>

namespace gsea {

    static double calculate_mean(const ExpressionData& expression,
                                 size_t gene_idx,
                                 const std::vector<size_t>& sample_indices) {
        double sum = 0.0;
        for (size_t col : sample_indices) {
            sum += expression.values()(gene_idx, col);
        }
        return sum / sample_indices.size();
    }

    std::vector<size_t> compute_gene_rank(const ExpressionData& expression,
                                           const std::vector<size_t>& disease_indices,
                                           const std::vector<size_t>& healthy_indices) {
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
            double diff = disease_mean - healthy_mean;
            gene_diffs.emplace_back(gene_idx, diff);
        }

        // Sort by difference in descending order
        std::sort(gene_diffs.begin(), gene_diffs.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        // Extract indices
        std::vector<size_t> ranked_indices;
        ranked_indices.reserve(num_genes);
        for (const auto& [idx, _] : gene_diffs) {
            ranked_indices.push_back(idx);
        }

        return ranked_indices;
    }

} // namespace gsea
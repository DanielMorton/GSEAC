#include "gsea/statistics.h"
#include "gsea/ranking.h"
#include "gsea/enrichment.h"
#include <random>
#include <algorithm>
#include <stdexcept>

#ifdef USE_PARALLEL_STL
#include <execution>
#endif

namespace gsea {

std::vector<size_t> generate_random_gene_rank(const ExpressionData& expression,
                                               size_t disease_size) {
    size_t num_cols = expression.num_samples();

    if (disease_size >= num_cols) {
        throw std::invalid_argument("Disease size must be less than total number of samples");
    }

    // Create and shuffle indices
    std::vector<size_t> all_indices(num_cols);
    std::iota(all_indices.begin(), all_indices.end(), size_t{0});

    std::random_device rd;
    std::mt19937 gen(rd());
    std::ranges::shuffle(all_indices, gen);

    // Split into disease and healthy
    auto disease_indices = std::span{all_indices.data(), disease_size};
    auto healthy_indices = std::span{all_indices.data() + disease_size, num_cols - disease_size};

    return compute_gene_rank(expression, disease_indices, healthy_indices);
}

std::vector<std::vector<double>> compute_null_distribution(
    const ExpressionData& expression,
    std::span<const GeneSet> gene_sets,
    size_t disease_size,
    size_t sample_size) {

    std::vector<size_t> indices(sample_size);
    std::iota(indices.begin(), indices.end(), size_t{0});

    std::vector<std::vector<double>> distribution(sample_size);

#ifdef USE_PARALLEL_STL
    std::transform(std::execution::par_unseq,
                   indices.begin(), indices.end(),
                   distribution.begin(),
                   [&](size_t) {
#else
    std::transform(indices.begin(), indices.end(), distribution.begin(), [&](size_t) {
#endif
        auto random_rank = generate_random_gene_rank(expression, disease_size);

        std::vector<double> scores;
        scores.reserve(gene_sets.size());
        for (const auto& gene_set : gene_sets) {
            scores.push_back(calculate_enrichment_score(gene_set, random_rank));
        }
        return scores;
    });

    return distribution;
}

std::vector<size_t> find_significant_sets(
    std::span<const double> actual_scores,
    std::span<const std::vector<double>> null_distribution,
    double p_value,
    size_t num_sets) {

    size_t sample_size = null_distribution.size();
    double corrected_p = p_value / num_sets;

    std::vector<size_t> significant;

    for (size_t i = 0; i < num_sets; ++i) {
        double actual_score = actual_scores[i];

        size_t count_greater = std::ranges::count_if(null_distribution,
            [&](const auto& sample) { return sample[i] >= actual_score; });

        double empirical_p = static_cast<double>(count_greater) / sample_size;
        if (empirical_p < corrected_p) {
            significant.push_back(i);
        }
    }

    return significant;
}

} // namespace gsea
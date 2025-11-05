#include "gsea/statistics.h"
#include "gsea/ranking.h"
#include "gsea/enrichment.h"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <ranges>

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
    auto all_indices = std::views::iota(size_t{0}, num_cols) | std::ranges::to<std::vector>();

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

    auto indices = std::views::iota(size_t{0}, sample_size) | std::ranges::to<std::vector>();

    std::vector<std::vector<double>> distribution(sample_size);

#ifdef USE_PARALLEL_STL
    std::transform(std::execution::par_unseq,
                   indices.begin(), indices.end(),
                   distribution.begin(),
                   [&](size_t) {
#else
    std::ranges::transform(indices, distribution.begin(), [&](size_t) {
#endif
        auto random_rank = generate_random_gene_rank(expression, disease_size);

        return gene_sets
            | std::views::transform([&](const auto& gene_set) {
                return calculate_enrichment_score(gene_set, random_rank);
            })
            | std::ranges::to<std::vector>();
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

    return std::views::iota(size_t{0}, num_sets)
        | std::views::filter([&](size_t i) {
            double actual_score = actual_scores[i];

            size_t count_greater = std::ranges::count_if(null_distribution,
                [&](const auto& sample) { return sample[i] >= actual_score; });

            double empirical_p = static_cast<double>(count_greater) / sample_size;
            return empirical_p < corrected_p;
        })
        | std::ranges::to<std::vector>();
}

} // namespace gsea
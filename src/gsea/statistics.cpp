#include "gsea/statistics.h"
#include "gsea/ranking.h"
#include "gsea/enrichment.h"
#include <random>
#include <algorithm>
#include <stdexcept>

namespace gsea {

std::vector<size_t> generate_random_gene_rank(const ExpressionData& expression,
                                               size_t disease_size) {
    size_t num_cols = expression.num_samples();
    
    if (disease_size >= num_cols) {
        throw std::invalid_argument("Disease size must be less than total number of samples");
    }

    // Create and shuffle indices
    std::vector<size_t> all_indices(num_cols);
    for (size_t i = 0; i < num_cols; ++i) {
        all_indices[i] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(all_indices.begin(), all_indices.end(), gen);

    // Split into disease and healthy
    std::vector<size_t> disease_indices(all_indices.begin(), 
                                        all_indices.begin() + disease_size);
    std::vector<size_t> healthy_indices(all_indices.begin() + disease_size, 
                                        all_indices.end());

    return compute_gene_rank(expression, disease_indices, healthy_indices);
}

std::vector<std::vector<double>> compute_null_distribution(
    const ExpressionData& expression,
    const std::vector<GeneSet>& gene_sets,
    size_t disease_size,
    size_t sample_size) {
    
    std::vector<std::vector<double>> distribution;
    distribution.reserve(sample_size);

    for (size_t i = 0; i < sample_size; ++i) {
        auto random_rank = generate_random_gene_rank(expression, disease_size);
        
        std::vector<double> scores;
        scores.reserve(gene_sets.size());
        
        for (const auto& gene_set : gene_sets) {
            scores.push_back(calculate_enrichment_score(gene_set, random_rank));
        }
        
        distribution.push_back(std::move(scores));
    }

    return distribution;
}

std::vector<size_t> find_significant_sets(
    const std::vector<double>& actual_scores,
    const std::vector<std::vector<double>>& null_distribution,
    double p_value,
    size_t num_sets) {
    
    size_t sample_size = null_distribution.size();
    double corrected_p = p_value / num_sets;

    std::vector<size_t> significant;

    for (size_t i = 0; i < num_sets; ++i) {
        double actual_score = actual_scores[i];
        
        size_t count_greater = 0;
        for (const auto& sample : null_distribution) {
            if (sample[i] >= actual_score) {
                ++count_greater;
            }
        }

        double empirical_p = static_cast<double>(count_greater) / sample_size;
        
        if (empirical_p < corrected_p) {
            significant.push_back(i);
        }
    }

    return significant;
}

} // namespace gsea
#include "gsea/analyzer.h"
#include "data_loader/expression_loader.h"
#include "data_loader/sample_loader.h"
#include "data_loader/geneset_loader.h"
#include "gsea/ranking.h"
#include "gsea/enrichment.h"
#include "gsea/statistics.h"
#include <iostream>
#include <stdexcept>
#include <format>
#include <ranges>

namespace gsea {

void GSEAAnalyzer::load_data(const std::string& exp_file,
                              const std::string& samp_file,
                              const std::string& geneset_file) {
    std::cout << "  Loading expression data...\n";
    expression_ = load_expression_data(exp_file);
    std::cout << std::format("    Loaded {} genes across {} samples\n",
              expression_->num_genes(), expression_->num_samples());

    std::cout << "  Loading sample data...\n";
    samples_ = load_sample_data(samp_file);
    std::cout << std::format("    Loaded {} samples ({} diseased, {} healthy)\n",
              samples_->num_samples(), samples_->num_diseased(), samples_->num_healthy());

    // Build sample name to column index mapping using ranges
    sample_to_column_.clear();
    for (auto [i, name] : std::views::enumerate(expression_->sample_names())) {
        sample_to_column_[std::string(name)] = i;
    }

    std::cout << "  Loading gene sets...\n";
    auto gene_names_vec = expression_->gene_names() | std::ranges::to<std::vector>();
    gene_sets_ = load_gene_sets(geneset_file, gene_names_vec);
    std::cout << std::format("    Loaded {} gene sets\n", gene_sets_.size());
}

std::vector<std::string> GSEAAnalyzer::get_gene_rank_order() {
    if (!samples_ || !expression_) {
        throw std::runtime_error("Data not loaded");
    }

    auto [disease_cols, healthy_cols] = [this]() {
        std::vector<size_t> disease, healthy;

        for (auto [i, status] : std::views::enumerate(samples_->disease_status())) {
            if (auto it = sample_to_column_.find(std::string(samples_->sample_names()[i]));
                it != sample_to_column_.end()) {
                (status == 1 ? disease : healthy).push_back(it->second);
            }
        }

        return std::pair{std::move(disease), std::move(healthy)};
    }();

    if (disease_cols.empty() || healthy_cols.empty()) {
        throw std::runtime_error(
            "No matching samples found between sample and expression files");
    }

    gene_rank_ = compute_gene_rank(*expression_, disease_cols, healthy_cols);

    return *gene_rank_
        | std::views::transform([this](size_t idx) {
            return std::string(expression_->gene_names()[idx]);
        })
        | std::ranges::to<std::vector>();
}

double GSEAAnalyzer::get_enrichment_score(const GeneSet& gene_set,
                                          const std::vector<size_t>* gene_rank) {
    if (!gene_rank) {
        if (!gene_rank_) {
            get_gene_rank_order();
        }
        gene_rank = &(*gene_rank_);
    }

    return calculate_enrichment_score(gene_set, *gene_rank);
}

std::unordered_map<std::string, double> GSEAAnalyzer::compute_all_enrichment_scores() {
    if (!gene_rank_) {
        get_gene_rank_order();
    }

    return gene_sets_
        | std::views::transform([this](const auto& gene_set) {
            return std::pair{
                std::string(gene_set.get_name()),
                calculate_enrichment_score(gene_set, *gene_rank_)
            };
        })
        | std::ranges::to<std::unordered_map>();
}

std::vector<std::string> GSEAAnalyzer::get_significant_sets(double p_value,
                                                             size_t sample_size) {
    if (!gene_rank_) {
        get_gene_rank_order();
    }

    if (!expression_ || !samples_) {
        throw std::runtime_error("Data not loaded");
    }

    // Compute actual enrichment scores
    auto actual_scores = gene_sets_
        | std::views::transform([this](const auto& gene_set) {
            return calculate_enrichment_score(gene_set, *gene_rank_);
        })
        | std::ranges::to<std::vector>();

    std::cout << std::format("  Generating null distribution with {} permutations...\n", sample_size);

    // Compute null distribution
    auto null_distribution = compute_null_distribution(
        *expression_,
        gene_sets_,
        samples_->num_diseased(),
        sample_size
    );

    // Find significant sets
    auto significant_indices = find_significant_sets(
        actual_scores,
        null_distribution,
        p_value,
        gene_sets_.size()
    );

    return significant_indices
        | std::views::transform([this](size_t idx) {
            return std::string(gene_sets_[idx].get_name());
        })
        | std::ranges::to<std::vector>();
}

bool GSEAAnalyzer::is_loaded() const {
    return expression_.has_value() && 
           samples_.has_value() && 
           !gene_sets_.empty();
}

} // namespace gsea
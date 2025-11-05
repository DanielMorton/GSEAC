#include "gsea/analyzer.h"
#include "data_loader/expression_loader.h"
#include "data_loader/sample_loader.h"
#include "data_loader/geneset_loader.h"
#include "gsea/ranking.h"
#include "gsea/enrichment.h"
#include "gsea/statistics.h"
#include <iostream>
#include <stdexcept>

namespace gsea {

void GSEAAnalyzer::load_data(const std::string& exp_file,
                              const std::string& samp_file,
                              const std::string& geneset_file) {
    std::cout << "  Loading expression data...\n";
    expression_ = load_expression_data(exp_file);
    std::cout << "    Loaded " << expression_->num_genes() 
              << " genes across " << expression_->num_samples() << " samples\n";

    std::cout << "  Loading sample data...\n";
    samples_ = load_sample_data(samp_file);
    std::cout << "    Loaded " << samples_->num_samples() << " samples ("
              << samples_->num_diseased() << " diseased, "
              << samples_->num_healthy() << " healthy)\n";

    // Build sample name to column index mapping
    sample_to_column_.clear();
    for (size_t i = 0; i < expression_->sample_names().size(); ++i) {
        sample_to_column_[expression_->sample_names()[i]] = i;
    }

    std::cout << "  Loading gene sets...\n";
    gene_sets_ = load_gene_sets(geneset_file, expression_->gene_names());
    std::cout << "    Loaded " << gene_sets_.size() << " gene sets\n";
}

std::vector<std::string> GSEAAnalyzer::get_gene_rank_order() {
    if (!samples_ || !expression_) {
        throw std::runtime_error("Data not loaded");
    }

    std::vector<size_t> disease_cols;
    std::vector<size_t> healthy_cols;

    for (size_t i = 0; i < samples_->disease_status().size(); ++i) {
        const auto& sample_name = samples_->sample_names()[i];
        auto it = sample_to_column_.find(sample_name);
        
        if (it != sample_to_column_.end()) {
            if (samples_->disease_status()[i] == 1) {
                disease_cols.push_back(it->second);
            } else {
                healthy_cols.push_back(it->second);
            }
        }
    }

    if (disease_cols.empty() || healthy_cols.empty()) {
        throw std::runtime_error(
            "No matching samples found between sample and expression files");
    }

    gene_rank_ = compute_gene_rank(*expression_, disease_cols, healthy_cols);

    std::vector<std::string> gene_names;
    gene_names.reserve(gene_rank_->size());
    for (size_t idx : *gene_rank_) {
        gene_names.push_back(expression_->gene_names()[idx]);
    }

    return gene_names;
}

double GSEAAnalyzer::get_enrichment_score(const GeneSet& gene_set,
                                          const std::vector<size_t>* gene_rank) {
    const std::vector<size_t>* rank = gene_rank;
    
    if (!rank) {
        if (!gene_rank_) {
            get_gene_rank_order();
        }
        rank = &(*gene_rank_);
    }

    return calculate_enrichment_score(gene_set, *rank);
}

std::unordered_map<std::string, double> GSEAAnalyzer::compute_all_enrichment_scores() {
    if (!gene_rank_) {
        get_gene_rank_order();
    }

    std::unordered_map<std::string, double> scores;
    for (const auto& gene_set : gene_sets_) {
        double score = calculate_enrichment_score(gene_set, *gene_rank_);
        scores[gene_set.get_name()] = score;
    }

    return scores;
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
    std::vector<double> actual_scores;
    actual_scores.reserve(gene_sets_.size());
    for (const auto& gene_set : gene_sets_) {
        actual_scores.push_back(calculate_enrichment_score(gene_set, *gene_rank_));
    }

    std::cout << "  Generating null distribution with " << sample_size 
              << " permutations...\n";

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

    std::vector<std::string> significant_names;
    significant_names.reserve(significant_indices.size());
    for (size_t idx : significant_indices) {
        significant_names.push_back(gene_sets_[idx].get_name());
    }

    return significant_names;
}

bool GSEAAnalyzer::is_loaded() const {
    return expression_.has_value() && 
           samples_.has_value() && 
           !gene_sets_.empty();
}

} // namespace gsea
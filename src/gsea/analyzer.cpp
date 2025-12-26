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

using namespace std;

namespace gsea {

GSEAAnalyzer::GSEAAnalyzer(const string& exp_file,
                           const string& samp_file,
                           const string& geneset_file)
    : expression_(load_expression_data(exp_file)),
      samples_(load_sample_data(samp_file))
{
    cout << "  Loading expression data...\n";
    cout << format("    Loaded {} genes across {} samples\n",
              expression_.num_genes(), expression_.num_samples());

    cout << "  Loading sample data...\n";
    cout << format("    Loaded {} samples ({} diseased, {} healthy)\n",
              samples_.num_samples(), samples_.num_diseased(), samples_.num_healthy());

    // Build sample name to column index mapping
    for (size_t i = 0; i < expression_.sample_names().size(); ++i) {
        sample_to_column_[string(expression_.sample_names()[i])] = i;
    }

    cout << "  Loading gene sets...\n";
    vector<string> gene_names_vec(expression_.gene_names().begin(),
                                             expression_.gene_names().end());
    gene_sets_ = load_gene_sets(geneset_file, gene_names_vec);
    cout << format("    Loaded {} gene sets\n", gene_sets_.size());
}

vector<string> GSEAAnalyzer::get_gene_rank_order() {
    vector<size_t> disease_cols;
    vector<size_t> healthy_cols;

    for (size_t i = 0; i < samples_.disease_status().size(); ++i) {
        auto sample_name = string(samples_.sample_names()[i]);
        if (auto it = sample_to_column_.find(sample_name); it != sample_to_column_.end()) {
            if (samples_.disease_status()[i] == 1) {
                disease_cols.push_back(it->second);
            } else {
                healthy_cols.push_back(it->second);
            }
        }
    }

    if (disease_cols.empty() || healthy_cols.empty()) {
        throw runtime_error(
            "No matching samples found between sample and expression files");
    }

    gene_rank_ = compute_gene_rank(expression_, disease_cols, healthy_cols);

    vector<string> gene_names;
    gene_names.reserve(gene_rank_.size());
    for (size_t idx : gene_rank_) {
        gene_names.emplace_back(expression_.gene_names()[idx]);
    }

    return gene_names;
}

double GSEAAnalyzer::get_enrichment_score(const GeneSet& gene_set,
                                          const vector<size_t>* gene_rank) {
    if (!gene_rank) {
        if (gene_rank_.empty()) {
            get_gene_rank_order();
        }
        gene_rank = &gene_rank_;
    }

    return calculate_enrichment_score(gene_set, *gene_rank);
}

unordered_map<string, double> GSEAAnalyzer::compute_all_enrichment_scores() {
    if (gene_rank_.empty()) {
        get_gene_rank_order();
    }

    unordered_map<string, double> scores;
    for (const auto& gene_set : gene_sets_) {
        scores[string(gene_set.get_name())] = calculate_enrichment_score(gene_set, gene_rank_);
    }

    return scores;
}

vector<string> GSEAAnalyzer::get_significant_sets(double p_value,
                                                             size_t sample_size) {
    if (gene_rank_.empty()) {
        get_gene_rank_order();
    }

    // Compute actual enrichment scores
    vector<double> actual_scores;
    actual_scores.reserve(gene_sets_.size());
    for (const auto& gene_set : gene_sets_) {
        actual_scores.push_back(calculate_enrichment_score(gene_set, gene_rank_));
    }

    cout << format("  Generating null distribution with {} permutations...\n", sample_size);

    // Compute null distribution
    auto null_distribution = compute_null_distribution(
        expression_,
        gene_sets_,
        samples_.num_diseased(),
        sample_size
    );

    // Find significant sets
    auto significant_indices = find_significant_sets(
        actual_scores,
        null_distribution,
        p_value,
        gene_sets_.size()
    );

    vector<string> significant_names;
    significant_names.reserve(significant_indices.size());
    for (size_t idx : significant_indices) {
        significant_names.emplace_back(gene_sets_[idx].get_name());
    }

    return significant_names;
}

} // namespace gsea
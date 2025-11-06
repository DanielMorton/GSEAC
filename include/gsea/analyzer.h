#pragma once

#include "types/expression_data.h"
#include "types/sample_data.h"
#include "types/gene_set.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <optional>

using namespace std;

namespace gsea {

class GSEAAnalyzer {
public:
    GSEAAnalyzer() = default;

    void load_data(const string& exp_file,
                   const string& samp_file,
                   const string& geneset_file);

    vector<string> get_gene_rank_order();

    double get_enrichment_score(const GeneSet& gene_set,
                               const vector<size_t>* gene_rank = nullptr);

    unordered_map<string, double> compute_all_enrichment_scores();

    vector<string> get_significant_sets(double p_value, size_t sample_size);

    size_t num_gene_sets() const { return gene_sets_.size(); }
    bool is_loaded() const;

private:
    optional<ExpressionData> expression_;
    optional<SampleData> samples_;
    vector<GeneSet> gene_sets_;
    optional<vector<size_t>> gene_rank_;
    unordered_map<string, size_t> sample_to_column_;
};

} // namespace gsea
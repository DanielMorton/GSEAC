#pragma once

#include "types/expression_data.h"
#include "types/sample_data.h"
#include "types/gene_set.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <optional>

namespace gsea {

    class GSEAAnalyzer {
    public:
        GSEAAnalyzer() = default;

        void load_data(const std::string& exp_file,
                       const std::string& samp_file,
                       const std::string& geneset_file);

        std::vector<std::string> get_gene_rank_order();

        double get_enrichment_score(const GeneSet& gene_set,
                                   const std::vector<size_t>* gene_rank = nullptr);

        std::unordered_map<std::string, double> compute_all_enrichment_scores();

        std::vector<std::string> get_significant_sets(double p_value, size_t sample_size);

        size_t num_gene_sets() const { return gene_sets_.size(); }
        bool is_loaded() const;

    private:
        std::optional<ExpressionData> expression_;
        std::optional<SampleData> samples_;
        std::vector<GeneSet> gene_sets_;
        std::optional<std::vector<size_t>> gene_rank_;
        std::unordered_map<std::string, size_t> sample_to_column_;
    };

} // namespace gsea
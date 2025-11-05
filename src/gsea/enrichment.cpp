#include "gsea/enrichment.h"
#include <limits>

namespace gsea {

    Eigen::VectorXd compute_brownian_bridge(const GeneSet& gene_set,
                                             const std::vector<size_t>& gene_rank) {
        size_t num_genes = gene_rank.size();
        Eigen::VectorXd bridge(num_genes);
        double cumsum = 0.0;

        for (size_t i = 0; i < num_genes; ++i) {
            size_t gene_idx = gene_rank[i];
            cumsum += gene_set.scores()(gene_idx);
            bridge(i) = cumsum;
        }

        return bridge;
    }

    double calculate_enrichment_score(const GeneSet& gene_set,
                                       const std::vector<size_t>& gene_rank) {
        auto bridge = compute_brownian_bridge(gene_set, gene_rank);

        double max_val = std::numeric_limits<double>::lowest();
        for (Eigen::Index i = 0; i < bridge.size(); ++i) {
            if (bridge(i) > max_val) {
                max_val = bridge(i);
            }
        }

        return max_val;
    }

} // namespace gsea
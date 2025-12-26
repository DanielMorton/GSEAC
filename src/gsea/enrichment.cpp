#include "gsea/enrichment.h"
#include <ranges>
#include <algorithm>

using namespace std;

namespace gsea {

Eigen::VectorXd compute_brownian_bridge(const GeneSet& gene_set,
                                         span<const size_t> gene_rank) {
    size_t num_genes = gene_rank.size();
    Eigen::VectorXd bridge(num_genes);

    double cumsum = 0.0;
    for (size_t i = 0; i < num_genes; ++i) {
        cumsum += gene_set.scores()(gene_rank[i]);
        bridge(i) = cumsum;
    }

    return bridge;
}

double calculate_enrichment_score(const GeneSet& gene_set,
                                   span<const size_t> gene_rank) {
    auto bridge = compute_brownian_bridge(gene_set, gene_rank);

    return ranges::max(bridge);
}

} // namespace gsea
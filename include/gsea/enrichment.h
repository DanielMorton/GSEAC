#pragma once

#include "types/gene_set.h"
#include <Eigen/Dense>
#include <span>

using namespace std;

namespace gsea {

    [[nodiscard]] Eigen::VectorXd compute_brownian_bridge(
        const GeneSet& gene_set,
        span<const size_t> gene_rank);

    [[nodiscard]] double calculate_enrichment_score(
        const GeneSet& gene_set,
        span<const size_t> gene_rank);

} // namespace gsea
#pragma once

#include "types/gene_set.h"
#include <Eigen/Dense>
#include <vector>

namespace gsea {

    Eigen::VectorXd compute_brownian_bridge(const GeneSet& gene_set,
                                             const std::vector<size_t>& gene_rank);

    double calculate_enrichment_score(const GeneSet& gene_set,
                                       const std::vector<size_t>& gene_rank);

} // namespace gsea
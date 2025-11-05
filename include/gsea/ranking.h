#pragma once

#include "types/expression_data.h"
#include <vector>

namespace gsea {

    std::vector<size_t> compute_gene_rank(const ExpressionData& expression,
                                           const std::vector<size_t>& disease_indices,
                                           const std::vector<size_t>& healthy_indices);

} // namespace gsea
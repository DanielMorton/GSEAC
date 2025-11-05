#pragma once

#include "types/expression_data.h"
#include <vector>
#include <span>

namespace gsea {

[[nodiscard]] std::vector<size_t> compute_gene_rank(
    const ExpressionData& expression,
    std::span<const size_t> disease_indices,
    std::span<const size_t> healthy_indices);

} // namespace gsea
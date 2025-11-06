#pragma once

#include "types/expression_data.h"
#include <vector>
#include <span>

using namespace std;

namespace gsea {

[[nodiscard]] vector<size_t> compute_gene_rank(
    const ExpressionData& expression,
    span<const size_t> disease_indices,
    span<const size_t> healthy_indices);

} // namespace gsea
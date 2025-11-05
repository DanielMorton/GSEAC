#include "types/gene_set.h"

namespace gsea {

GeneSet::GeneSet(std::string name, size_t gene_count, Eigen::VectorXd scores)
    : name_(std::move(name))
    , gene_count_(gene_count)
    , scores_(std::move(scores)) {}

std::optional<double> GeneSet::get_score(size_t gene_idx) const noexcept {
    if (gene_idx < static_cast<size_t>(scores_.size())) {
        return scores_(gene_idx);
    }
    return std::nullopt;
}

} // namespace gsea
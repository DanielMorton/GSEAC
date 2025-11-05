#include "types/gene_set.h"

using namespace std;

namespace gsea {

GeneSet::GeneSet(string name, size_t gene_count, Eigen::VectorXd scores)
    : name_(std::move(name))
    , gene_count_(gene_count)
    , scores_(std::move(scores)) {}

optional<double> GeneSet::get_score(size_t gene_idx) const noexcept {
    if (gene_idx < static_cast<size_t>(scores_.size())) {
        return scores_(gene_idx);
    }
    return nullopt;
}

} // namespace gsea
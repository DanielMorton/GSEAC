#include "types/expression_data.h"

using namespace std;

namespace gsea {

ExpressionData::ExpressionData(Eigen::MatrixXd values,
                               vector<string> gene_names,
                               vector<string> sample_names)
    : values_(std::move(values))
    , gene_names_(std::move(gene_names))
    , sample_names_(std::move(sample_names)) {}

optional<string_view> ExpressionData::gene_name(size_t index) const noexcept {
    if (index < gene_names_.size()) {
        return gene_names_[index];
    }
    return nullopt;
}

optional<string_view> ExpressionData::sample_name(size_t index) const noexcept {
    if (index < sample_names_.size()) {
        return sample_names_[index];
    }
    return nullopt;
}

optional<double> ExpressionData::get_value(size_t gene_idx, size_t sample_idx) const noexcept {
    if (gene_idx < num_genes() && sample_idx < num_samples()) {
        return values_(gene_idx, sample_idx);
    }
    return nullopt;
}

} // namespace gsea
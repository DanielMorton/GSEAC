#include "types/expression_data.h"

namespace gsea {

ExpressionData::ExpressionData(Eigen::MatrixXd values,
                               std::vector<std::string> gene_names,
                               std::vector<std::string> sample_names)
    : values_(std::move(values))
    , gene_names_(std::move(gene_names))
    , sample_names_(std::move(sample_names)) {}

std::optional<std::string_view> ExpressionData::gene_name(size_t index) const noexcept {
    if (index < gene_names_.size()) {
        return gene_names_[index];
    }
    return std::nullopt;
}

std::optional<std::string_view> ExpressionData::sample_name(size_t index) const noexcept {
    if (index < sample_names_.size()) {
        return sample_names_[index];
    }
    return std::nullopt;
}

std::optional<double> ExpressionData::get_value(size_t gene_idx, size_t sample_idx) const noexcept {
    if (gene_idx < num_genes() && sample_idx < num_samples()) {
        return values_(gene_idx, sample_idx);
    }
    return std::nullopt;
}

} // namespace gsea
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <optional>
#include <span>

namespace gsea {

class ExpressionData {
public:
    ExpressionData(Eigen::MatrixXd values,
                   std::vector<std::string> gene_names,
                   std::vector<std::string> sample_names);

    [[nodiscard]] size_t num_genes() const noexcept { return gene_names_.size(); }
    [[nodiscard]] size_t num_samples() const noexcept { return sample_names_.size(); }

    [[nodiscard]] std::optional<std::string_view> gene_name(size_t index) const noexcept;
    [[nodiscard]] std::optional<std::string_view> sample_name(size_t index) const noexcept;
    [[nodiscard]] std::optional<double> get_value(size_t gene_idx, size_t sample_idx) const noexcept;

    [[nodiscard]] const Eigen::MatrixXd& values() const noexcept { return values_; }
    [[nodiscard]] std::span<const std::string> gene_names() const noexcept { return gene_names_; }
    [[nodiscard]] std::span<const std::string> sample_names() const noexcept { return sample_names_; }

private:
    Eigen::MatrixXd values_;
    std::vector<std::string> gene_names_;
    std::vector<std::string> sample_names_;
};

} // namespace gsea
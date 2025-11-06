#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <optional>
#include <span>

using namespace std;

namespace gsea {

class ExpressionData {
public:
    ExpressionData(Eigen::MatrixXd values,
                   vector<string> gene_names,
                   vector<string> sample_names);

    [[nodiscard]] size_t num_genes() const noexcept { return gene_names_.size(); }
    [[nodiscard]] size_t num_samples() const noexcept { return sample_names_.size(); }

    [[nodiscard]] optional<string_view> gene_name(size_t index) const noexcept;
    [[nodiscard]] optional<string_view> sample_name(size_t index) const noexcept;
    [[nodiscard]] optional<double> get_value(size_t gene_idx, size_t sample_idx) const noexcept;

    [[nodiscard]] const Eigen::MatrixXd& values() const noexcept { return values_; }
    [[nodiscard]] span<const string> gene_names() const noexcept { return gene_names_; }
    [[nodiscard]] span<const string> sample_names() const noexcept { return sample_names_; }

private:
    Eigen::MatrixXd values_;
    vector<string> gene_names_;
    vector<string> sample_names_;
};

} // namespace gsea
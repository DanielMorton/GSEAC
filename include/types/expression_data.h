#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <optional>

namespace gsea {

    class ExpressionData {
    public:
        ExpressionData(Eigen::MatrixXd values,
                       std::vector<std::string> gene_names,
                       std::vector<std::string> sample_names);

        size_t num_genes() const { return gene_names_.size(); }
        size_t num_samples() const { return sample_names_.size(); }

        std::optional<std::string> gene_name(size_t index) const;
        std::optional<std::string> sample_name(size_t index) const;
        std::optional<double> get_value(size_t gene_idx, size_t sample_idx) const;

        const Eigen::MatrixXd& values() const { return values_; }
        const std::vector<std::string>& gene_names() const { return gene_names_; }
        const std::vector<std::string>& sample_names() const { return sample_names_; }

    private:
        Eigen::MatrixXd values_;
        std::vector<std::string> gene_names_;
        std::vector<std::string> sample_names_;
    };

} // namespace gsea
#pragma once

#include <Eigen/Dense>
#include <string>
#include <optional>

namespace gsea {

    class GeneSet {
    public:
        GeneSet(std::string name, size_t gene_count, Eigen::VectorXd scores);

        size_t size() const { return gene_count_; }
        std::optional<double> get_score(size_t gene_idx) const;
        const std::string& get_name() const { return name_; }

        const Eigen::VectorXd& scores() const { return scores_; }

    private:
        std::string name_;
        size_t gene_count_;
        Eigen::VectorXd scores_;
    };

} // namespace gsea
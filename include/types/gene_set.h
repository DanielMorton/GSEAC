#pragma once

#include <Eigen/Dense>
#include <string>
#include <optional>
#include <string_view>

using namespace std;

namespace gsea {

class GeneSet {
public:
    GeneSet(string name, size_t gene_count, Eigen::VectorXd scores);

    [[nodiscard]] size_t size() const noexcept { return gene_count_; }
    [[nodiscard]] optional<double> get_score(size_t gene_idx) const noexcept;
    [[nodiscard]] string_view get_name() const noexcept { return name_; }

    [[nodiscard]] const Eigen::VectorXd& scores() const noexcept { return scores_; }

private:
    string name_;
    size_t gene_count_;
    Eigen::VectorXd scores_;
};

} // namespace gsea
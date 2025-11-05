#include "data_loader/geneset_loader.h"
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <cmath>
#include <ranges>
#include <string_view>
#include <format>

namespace gsea {

static auto split(std::string_view str, char delimiter) {
    return str
        | std::views::split(delimiter)
        | std::views::transform([](auto&& rng) {
            return std::string(rng.begin(), rng.end());
        })
        | std::ranges::to<std::vector>();
}

static std::string trim(std::string_view str) {
    constexpr std::string_view whitespace = " \t\r\n";
    auto start = str.find_first_not_of(whitespace);
    if (start == std::string_view::npos) return "";
    auto end = str.find_last_not_of(whitespace);
    return std::string(str.substr(start, end - start + 1));
}

std::vector<GeneSet> load_gene_sets(const std::string& filepath,
                                     const std::vector<std::string>& gene_names) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error(std::format("Failed to open gene set file: {}", filepath));
    }

    size_t num_genes = gene_names.size();
    std::vector<GeneSet> gene_sets;

    std::string line;
    size_t line_num = 0;

    while (std::getline(file, line)) {
        ++line_num;
        if (line.empty()) continue;

        auto tokens = split(line, '\t');
        if (tokens.size() < 3) {
            std::cerr << std::format("Warning: Skipping line {} in gene set file: insufficient columns\n",
                line_num);
            continue;
        }

        std::string set_name = trim(tokens[0]);

        // Build set of genes (skip name and description columns)
        auto genes_in_set = tokens
            | std::views::drop(2)
            | std::views::transform([](const auto& s) { return trim(s); })
            | std::views::filter([](const auto& s) { return !s.empty(); })
            | std::ranges::to<std::unordered_set<std::string>>();

        if (genes_in_set.empty()) {
            std::cerr << std::format("Warning: Skipping gene set '{}': contains no genes\n", set_name);
            continue;
        }

        // Create boolean mask
        auto gene_mask = gene_names
            | std::views::transform([&](const auto& gene) {
                return genes_in_set.contains(gene);
            })
            | std::ranges::to<std::vector>();

        size_t gene_count = std::ranges::count(gene_mask, true);

        if (gene_count == 0) {
            std::cerr << std::format("Warning: Skipping gene set '{}': no genes match expression data\n",
                set_name);
            continue;
        }

        // Precompute scores
        double up_score = std::sqrt(static_cast<double>(num_genes - gene_count) / gene_count);
        double down_score = -1.0 / up_score;

        auto scores = gene_mask
            | std::views::transform([&](bool in_set) {
                return in_set ? up_score : down_score;
            })
            | std::ranges::to<std::vector>();

        Eigen::VectorXd scores_vec = Eigen::Map<Eigen::VectorXd>(scores.data(), scores.size());

        gene_sets.emplace_back(std::move(set_name), gene_count, std::move(scores_vec));
    }

    if (gene_sets.empty()) {
        throw std::runtime_error("No valid gene sets found in file");
    }

    return gene_sets;
}

} // namespace gsea
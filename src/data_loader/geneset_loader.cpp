#include "data_loader/geneset_loader.h"
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <cmath>
#include <string_view>
#include <format>

namespace gsea {

static std::vector<std::string> split(std::string_view str, char delimiter) {
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = str.find(delimiter);

    while (end != std::string_view::npos) {
        tokens.emplace_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find(delimiter, start);
    }
    tokens.emplace_back(str.substr(start));

    return tokens;
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
        std::unordered_set<std::string> genes_in_set;
        for (size_t i = 2; i < tokens.size(); ++i) {
            auto gene = trim(tokens[i]);
            if (!gene.empty()) {
                genes_in_set.insert(std::move(gene));
            }
        }

        if (genes_in_set.empty()) {
            std::cerr << std::format("Warning: Skipping gene set '{}': contains no genes\n", set_name);
            continue;
        }

        // Create boolean mask and count genes
        std::vector<bool> gene_mask;
        gene_mask.reserve(num_genes);
        size_t gene_count = 0;

        for (const auto& gene : gene_names) {
            bool in_set = genes_in_set.contains(gene);
            gene_mask.push_back(in_set);
            if (in_set) ++gene_count;
        }

        if (gene_count == 0) {
            std::cerr << std::format("Warning: Skipping gene set '{}': no genes match expression data\n",
                set_name);
            continue;
        }

        // Precompute scores
        double up_score = std::sqrt(static_cast<double>(num_genes - gene_count) / gene_count);
        double down_score = -1.0 / up_score;

        Eigen::VectorXd scores(num_genes);
        for (size_t i = 0; i < num_genes; ++i) {
            scores(i) = gene_mask[i] ? up_score : down_score;
        }

        gene_sets.emplace_back(std::move(set_name), gene_count, std::move(scores));
    }

    if (gene_sets.empty()) {
        throw std::runtime_error("No valid gene sets found in file");
    }

    return gene_sets;
}

} // namespace gsea
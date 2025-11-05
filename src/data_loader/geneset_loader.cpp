#include "data_loader/geneset_loader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <iostream>
#include <cmath>

namespace gsea {

static std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

static std::string trim(const std::string& str) {
    auto start = str.find_first_not_of(" \t\r\n");
    auto end = str.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return str.substr(start, end - start + 1);
}

std::vector<GeneSet> load_gene_sets(const std::string& filepath,
                                     const std::vector<std::string>& gene_names) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Failed to open gene set file: " + filepath);
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
            std::cerr << "Warning: Skipping line " << line_num 
                      << " in gene set file: insufficient columns\n";
            continue;
        }

        std::string set_name = trim(tokens[0]);

        // Build set of genes (skip name and description columns)
        std::unordered_set<std::string> genes_in_set;
        for (size_t i = 2; i < tokens.size(); ++i) {
            std::string gene = trim(tokens[i]);
            if (!gene.empty()) {
                genes_in_set.insert(gene);
            }
        }

        if (genes_in_set.empty()) {
            std::cerr << "Warning: Skipping gene set '" << set_name 
                      << "': contains no genes\n";
            continue;
        }

        // Create boolean mask
        std::vector<bool> gene_mask;
        for (const auto& gene : gene_names) {
            gene_mask.push_back(genes_in_set.count(gene) > 0);
        }

        size_t gene_count = 0;
        for (bool in_set : gene_mask) {
            if (in_set) ++gene_count;
        }

        if (gene_count == 0) {
            std::cerr << "Warning: Skipping gene set '" << set_name 
                      << "': no genes match expression data\n";
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
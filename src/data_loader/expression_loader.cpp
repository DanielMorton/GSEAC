#include "data_loader/expression_loader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

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

ExpressionData load_expression_data(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Failed to open expression file: " + filepath);
    }

    std::string line;
    
    // Read header
    if (!std::getline(file, line)) {
        throw std::runtime_error("Expression file is empty");
    }

    auto headers = split(line, '\t');
    if (headers.size() < 2) {
        throw std::runtime_error("Expression file must have at least 2 columns");
    }

    // Extract sample names (skip first column "SYMBOL")
    std::vector<std::string> sample_names;
    for (size_t i = 1; i < headers.size(); ++i) {
        sample_names.push_back(trim(headers[i]));
    }
    size_t num_samples = sample_names.size();

    // Read data
    std::vector<std::string> gene_names;
    std::vector<double> values;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        auto tokens = split(line, '\t');
        if (tokens.size() != num_samples + 1) {
            throw std::runtime_error(
                "Inconsistent number of columns at gene: " + tokens[0]);
        }

        gene_names.push_back(trim(tokens[0]));

        for (size_t i = 1; i < tokens.size(); ++i) {
            try {
                values.push_back(std::stod(tokens[i]));
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "Failed to parse expression value: " + tokens[i]);
            }
        }
    }

    size_t num_genes = gene_names.size();
    if (num_genes == 0) {
        throw std::runtime_error("Expression file contains no gene rows");
    }

    // Create matrix
    Eigen::MatrixXd matrix(num_genes, num_samples);
    for (size_t i = 0; i < num_genes; ++i) {
        for (size_t j = 0; j < num_samples; ++j) {
            matrix(i, j) = values[i * num_samples + j];
        }
    }

    return ExpressionData(std::move(matrix), std::move(gene_names), std::move(sample_names));
}

} // namespace gsea
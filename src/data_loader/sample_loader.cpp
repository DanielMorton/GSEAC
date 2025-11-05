#include "data_loader/sample_loader.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

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

SampleData load_sample_data(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Failed to open sample file: " + filepath);
    }

    std::vector<std::string> sample_names;
    std::vector<uint8_t> disease_status;

    std::string line;
    size_t line_num = 0;

    while (std::getline(file, line)) {
        ++line_num;
        if (line.empty()) continue;

        auto tokens = split(line, '\t');
        if (tokens.size() != 2) {
            throw std::runtime_error(
                "Invalid format at line " + std::to_string(line_num) +
                ": expected 2 columns, found " + std::to_string(tokens.size()));
        }

        sample_names.push_back(trim(tokens[0]));

        try {
            int status = std::stoi(trim(tokens[1]));
            if (status != 0 && status != 1) {
                throw std::runtime_error(
                    "Invalid disease status at line " + std::to_string(line_num) +
                    ": must be 0 or 1, found " + std::to_string(status));
            }
            disease_status.push_back(static_cast<uint8_t>(status));
        } catch (const std::invalid_argument&) {
            throw std::runtime_error(
                "Invalid disease status at line " + std::to_string(line_num) +
                ": '" + tokens[1] + "' is not a valid number");
        }
    }

    if (sample_names.empty()) {
        throw std::runtime_error("Sample file contains no data");
    }

    return SampleData(std::move(sample_names), std::move(disease_status));
}

} // namespace gsea
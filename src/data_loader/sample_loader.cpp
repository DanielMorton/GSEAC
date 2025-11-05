#include "data_loader/sample_loader.h"
#include <fstream>
#include <stdexcept>
#include <string_view>
#include <format>

using namespace std;

namespace gsea {

static vector<string> split(string_view str, char delimiter) {
    vector<string> tokens;
    size_t start = 0;
    size_t end = str.find(delimiter);

    while (end != string_view::npos) {
        tokens.emplace_back(str.substr(start, end - start));
        start = end + 1;
        end = str.find(delimiter, start);
    }
    tokens.emplace_back(str.substr(start));

    return tokens;
}

static string trim(string_view str) {
    constexpr string_view whitespace = " \t\r\n";
    auto start = str.find_first_not_of(whitespace);
    if (start == string_view::npos) return "";
    auto end = str.find_last_not_of(whitespace);
    return string(str.substr(start, end - start + 1));
}

SampleData load_sample_data(const string& filepath) {
    ifstream file(filepath);
    if (!file) {
        throw runtime_error(format("Failed to open sample file: {}", filepath));
    }

    vector<string> sample_names;
    vector<uint8_t> disease_status;

    string line;
    size_t line_num = 0;

    while (getline(file, line)) {
        ++line_num;
        if (line.empty()) continue;

        auto tokens = split(line, '\t');
        if (tokens.size() != 2) {
            throw runtime_error(
                format("Invalid format at line {}: expected 2 columns, found {}",
                    line_num, tokens.size()));
        }

        sample_names.push_back(trim(tokens[0]));

        try {
            int status = stoi(trim(tokens[1]));
            if (status != 0 && status != 1) {
                throw runtime_error(
                    format("Invalid disease status at line {}: must be 0 or 1, found {}",
                        line_num, status));
            }
            disease_status.push_back(static_cast<uint8_t>(status));
        } catch (const invalid_argument&) {
            throw runtime_error(
                format("Invalid disease status at line {}: '{}' is not a valid number",
                    line_num, tokens[1]));
        }
    }

    if (sample_names.empty()) {
        throw runtime_error("Sample file contains no data");
    }

    return SampleData(std::move(sample_names), std::move(disease_status));
}

} // namespace gsea
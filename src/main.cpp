#include "gsea/analyzer.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <format>
#include <ranges>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << std::format("Usage: {} <expression_file> <sample_file> <geneset_file>\n", argv[0]);
        std::cerr << "Please specify an expression file, sample file, and gene set file.\n";
        return 1;
    }

    const auto [exp_file, samp_file, kegg_file] = std::tuple{argv[1], argv[2], argv[3]};

    try {
        std::cout << "Loading data...\n";
        gsea::GSEAAnalyzer analyzer;
        analyzer.load_data(exp_file, samp_file, kegg_file);

        std::cout << "Computing enrichment scores...\n";
        auto es_scores = analyzer.compute_all_enrichment_scores();

        // Sort by enrichment score (descending) using ranges
        std::vector<std::pair<std::string, double>> sorted_scores;
        sorted_scores.reserve(es_scores.size());
        std::ranges::copy(es_scores, std::back_inserter(sorted_scores));

        std::ranges::sort(sorted_scores, std::ranges::greater{}, &std::pair<std::string, double>::second);

        // Write to file
        std::ofstream file("kegg_enrichment_scores.txt");
        if (!file) {
            throw std::runtime_error("Failed to create output file");
        }

        for (const auto& [name, score] : sorted_scores) {
            file << std::format("{}\t{}\n", name, score);
        }

        std::cout << "Computing statistically significant gene sets...\n";
        auto sig_sets = analyzer.get_significant_sets(0.05, 100);

        std::cout << "Significant gene sets:\n";
        for (const auto& set_name : sig_sets) {
            std::cout << set_name << '\n';
        }

    } catch (const std::exception& e) {
        std::cerr << std::format("Error: {}\n", e.what());
        return 1;
    }

    return 0;
}
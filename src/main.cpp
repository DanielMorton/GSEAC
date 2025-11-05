#include "gsea/analyzer.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <expression_file> <sample_file> <geneset_file>\n";
        std::cerr << "Please specify an expression file, sample file, and gene set file.\n";
        return 1;
    }

    const char* exp_file = argv[1];
    const char* samp_file = argv[2];
    const char* kegg_file = argv[3];

    try {
        std::cout << "Loading data...\n";
        gsea::GSEAAnalyzer analyzer;
        analyzer.load_data(exp_file, samp_file, kegg_file);

        std::cout << "Computing enrichment scores...\n";
        auto es_scores = analyzer.compute_all_enrichment_scores();

        // Sort by enrichment score (descending)
        std::vector<std::pair<std::string, double>> sorted_scores(
            es_scores.begin(), es_scores.end()
        );
        std::sort(sorted_scores.begin(), sorted_scores.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        // Write to file
        std::ofstream file("kegg_enrichment_scores.txt");
        if (!file) {
            throw std::runtime_error("Failed to create output file");
        }
        
        for (const auto& [name, score] : sorted_scores) {
            file << name << "\t" << score << "\n";
        }
        file.close();

        std::cout << "Computing statistically significant gene sets...\n";
        auto sig_sets = analyzer.get_significant_sets(0.05, 100);

        std::cout << "Significant gene sets:\n";
        for (const auto& set_name : sig_sets) {
            std::cout << set_name << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
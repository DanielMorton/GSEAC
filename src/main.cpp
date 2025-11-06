#include "gsea/analyzer.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <format>
#include <ranges>

using namespace std;
using namespace gsea;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << format("Usage: {} <expression_file> <sample_file> <geneset_file>\n", argv[0]);
        cerr << "Please specify an expression file, sample file, and gene set file.\n";
        return 1;
    }

    const auto [exp_file, samp_file, kegg_file] = tuple{argv[1], argv[2], argv[3]};

    try {
        cout << "Loading data...\n";
        GSEAAnalyzer analyzer(exp_file, samp_file, kegg_file);

        cout << "Computing enrichment scores...\n";
        auto es_scores = analyzer.compute_all_enrichment_scores();

        // Sort by enrichment score (descending) using ranges
        vector<pair<string, double>> sorted_scores;
        sorted_scores.reserve(es_scores.size());
        ranges::copy(es_scores, back_inserter(sorted_scores));

        ranges::sort(sorted_scores, ranges::greater{}, &pair<string, double>::second);

        // Write to file
        ofstream file("kegg_enrichment_scores.txt");
        if (!file) {
            throw runtime_error("Failed to create output file");
        }

        for (const auto& [name, score] : sorted_scores) {
            file << format("{}\t{}\n", name, score);
        }

        cout << "Computing statistically significant gene sets...\n";
        auto sig_sets = analyzer.get_significant_sets(0.05, 100);

        cout << "Significant gene sets:\n";
        for (const auto& set_name : sig_sets) {
            cout << set_name << '\n';
        }

    } catch (const exception& e) {
        cerr << format("Error: {}\n", e.what());
        return 1;
    }

    return 0;
}
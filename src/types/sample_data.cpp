#include "types/sample_data.h"
#include <stdexcept>
#include <algorithm>
#include <ranges>

using namespace std;

namespace gsea {

SampleData::SampleData(vector<string> sample_names,
                       vector<uint8_t> disease_status)
    : sample_names_(std::move(sample_names))
    , disease_status_(std::move(disease_status)) {
    if (sample_names_.size() != disease_status_.size()) {
        throw invalid_argument(
            "Sample names and disease status must have the same length");
    }
}

size_t SampleData::num_diseased() const noexcept {
    return ranges::count(disease_status_, 1);
}

size_t SampleData::num_healthy() const noexcept {
    return ranges::count(disease_status_, 0);
}

vector<size_t> SampleData::get_disease_indices() const {
    vector<size_t> indices;
    for (size_t i = 0; i < disease_status_.size(); ++i) {
        if (disease_status_[i] == 1) {
            indices.push_back(i);
        }
    }
    return indices;
}

vector<size_t> SampleData::get_healthy_indices() const {
    vector<size_t> indices;
    for (size_t i = 0; i < disease_status_.size(); ++i) {
        if (disease_status_[i] == 0) {
            indices.push_back(i);
        }
    }
    return indices;
}

optional<bool> SampleData::is_diseased(size_t index) const noexcept {
    if (index < disease_status_.size()) {
        return disease_status_[index] == 1;
    }
    return nullopt;
}

optional<string_view> SampleData::get_sample_name(size_t index) const noexcept {
    if (index < sample_names_.size()) {
        return sample_names_[index];
    }
    return nullopt;
}

} // namespace gsea
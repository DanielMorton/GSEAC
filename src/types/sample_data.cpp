#include "types/sample_data.h"
#include <stdexcept>
#include <algorithm>

namespace gsea {

    SampleData::SampleData(std::vector<std::string> sample_names,
                           std::vector<uint8_t> disease_status)
        : sample_names_(std::move(sample_names))
        , disease_status_(std::move(disease_status)) {
        if (sample_names_.size() != disease_status_.size()) {
            throw std::invalid_argument(
                "Sample names and disease status must have the same length");
        }
    }

    size_t SampleData::num_diseased() const {
        return std::count(disease_status_.begin(), disease_status_.end(), 1);
    }

    size_t SampleData::num_healthy() const {
        return std::count(disease_status_.begin(), disease_status_.end(), 0);
    }

    std::vector<size_t> SampleData::get_disease_indices() const {
        std::vector<size_t> indices;
        for (size_t i = 0; i < disease_status_.size(); ++i) {
            if (disease_status_[i] == 1) {
                indices.push_back(i);
            }
        }
        return indices;
    }

    std::vector<size_t> SampleData::get_healthy_indices() const {
        std::vector<size_t> indices;
        for (size_t i = 0; i < disease_status_.size(); ++i) {
            if (disease_status_[i] == 0) {
                indices.push_back(i);
            }
        }
        return indices;
    }

    std::optional<bool> SampleData::is_diseased(size_t index) const {
        if (index < disease_status_.size()) {
            return disease_status_[index] == 1;
        }
        return std::nullopt;
    }

    std::optional<std::string> SampleData::get_sample_name(size_t index) const {
        if (index < sample_names_.size()) {
            return sample_names_[index];
        }
        return std::nullopt;
    }

} // namespace gsea
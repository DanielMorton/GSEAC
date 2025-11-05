#include "types/sample_data.h"
#include <stdexcept>
#include <algorithm>
#include <ranges>

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

size_t SampleData::num_diseased() const noexcept {
    return std::ranges::count(disease_status_, 1);
}

size_t SampleData::num_healthy() const noexcept {
    return std::ranges::count(disease_status_, 0);
}

std::vector<size_t> SampleData::get_disease_indices() const {
    auto indices = std::views::iota(size_t{0}, disease_status_.size())
                 | std::views::filter([this](size_t i) { return disease_status_[i] == 1; })
                 | std::ranges::to<std::vector>();
    return indices;
}

std::vector<size_t> SampleData::get_healthy_indices() const {
    auto indices = std::views::iota(size_t{0}, disease_status_.size())
                 | std::views::filter([this](size_t i) { return disease_status_[i] == 0; })
                 | std::ranges::to<std::vector>();
    return indices;
}

std::optional<bool> SampleData::is_diseased(size_t index) const noexcept {
    if (index < disease_status_.size()) {
        return disease_status_[index] == 1;
    }
    return std::nullopt;
}

std::optional<std::string_view> SampleData::get_sample_name(size_t index) const noexcept {
    if (index < sample_names_.size()) {
        return sample_names_[index];
    }
    return std::nullopt;
}

} // namespace gsea
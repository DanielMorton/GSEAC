#pragma once

#include <vector>
#include <string>
#include <optional>
#include <cstdint>
#include <span>

namespace gsea {

class SampleData {
public:
    SampleData(std::vector<std::string> sample_names,
               std::vector<uint8_t> disease_status);

    [[nodiscard]] size_t num_samples() const noexcept { return sample_names_.size(); }
    [[nodiscard]] size_t num_diseased() const noexcept;
    [[nodiscard]] size_t num_healthy() const noexcept;

    [[nodiscard]] std::vector<size_t> get_disease_indices() const;
    [[nodiscard]] std::vector<size_t> get_healthy_indices() const;

    [[nodiscard]] std::optional<bool> is_diseased(size_t index) const noexcept;
    [[nodiscard]] std::optional<std::string_view> get_sample_name(size_t index) const noexcept;

    [[nodiscard]] std::span<const std::string> sample_names() const noexcept { return sample_names_; }
    [[nodiscard]] std::span<const uint8_t> disease_status() const noexcept { return disease_status_; }

private:
    std::vector<std::string> sample_names_;
    std::vector<uint8_t> disease_status_;
};

} // namespace gsea
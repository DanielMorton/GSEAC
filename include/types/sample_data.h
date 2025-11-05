#pragma once

#include <vector>
#include <string>
#include <optional>
#include <cstdint>

namespace gsea {

    class SampleData {
    public:
        SampleData(std::vector<std::string> sample_names,
                   std::vector<uint8_t> disease_status);

        size_t num_samples() const { return sample_names_.size(); }
        size_t num_diseased() const;
        size_t num_healthy() const;

        std::vector<size_t> get_disease_indices() const;
        std::vector<size_t> get_healthy_indices() const;

        std::optional<bool> is_diseased(size_t index) const;
        std::optional<std::string> get_sample_name(size_t index) const;

        const std::vector<std::string>& sample_names() const { return sample_names_; }
        const std::vector<uint8_t>& disease_status() const { return disease_status_; }

    private:
        std::vector<std::string> sample_names_;
        std::vector<uint8_t> disease_status_;
    };

} // namespace gsea
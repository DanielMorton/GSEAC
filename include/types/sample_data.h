#pragma once

#include <vector>
#include <string>
#include <optional>
#include <span>

using namespace std;

namespace gsea {

class SampleData {
public:
    SampleData(vector<string> sample_names,
               vector<uint8_t> disease_status);

    [[nodiscard]] size_t num_samples() const noexcept { return sample_names_.size(); }
    [[nodiscard]] size_t num_diseased() const noexcept;
    [[nodiscard]] size_t num_healthy() const noexcept;

    [[nodiscard]] vector<size_t> get_disease_indices() const;
    [[nodiscard]] vector<size_t> get_healthy_indices() const;

    [[nodiscard]] optional<bool> is_diseased(size_t index) const noexcept;
    [[nodiscard]] optional<string_view> get_sample_name(size_t index) const noexcept;

    [[nodiscard]] span<const string> sample_names() const noexcept { return sample_names_; }
    [[nodiscard]] span<const uint8_t> disease_status() const noexcept { return disease_status_; }

private:
    vector<string> sample_names_;
    vector<uint8_t> disease_status_;
};

} // namespace gsea
#pragma once

#include "types/sample_data.h"
#include <string>

namespace gsea {

SampleData load_sample_data(const std::string& filepath);

} // namespace gsea
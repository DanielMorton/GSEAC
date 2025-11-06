#pragma once

#include "types/sample_data.h"
#include <string>

using namespace std;

namespace gsea {

SampleData load_sample_data(const string& filepath);

} // namespace gsea
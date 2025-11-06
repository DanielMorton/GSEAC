#pragma once

#include "types/expression_data.h"
#include <string>

using namespace std;

namespace gsea {

ExpressionData load_expression_data(const string& filepath);

} // namespace gsea
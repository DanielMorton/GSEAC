#pragma once

#include "types/expression_data.h"
#include <string>

namespace gsea {

    ExpressionData load_expression_data(const std::string& filepath);

} // namespace gsea
#pragma once

#include "types/gene_set.h"
#include <string>
#include <vector>

namespace gsea {

    std::vector<GeneSet> load_gene_sets(const std::string& filepath,
                                         const std::vector<std::string>& gene_names);

} // namespace gsea
#pragma once

#include "types/gene_set.h"
#include <string>
#include <vector>

using namespace std;

namespace gsea {

vector<GeneSet> load_gene_sets(const string& filepath,
                                     const vector<string>& gene_names);

} // namespace gsea
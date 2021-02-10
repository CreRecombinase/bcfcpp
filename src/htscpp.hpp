#pragma once

#include <vector>
#include <memory>
#include <cstddef>
#include "htslib/hts.h"
#include "htslib/bgzf.h"

int read_line_bcf(std::vector<std::byte> &out_vec, BGZF* fp);

#pragma once

#include <cstdint>
#include <vector>


class Chunk{
  std::uint64_t chunk_start; //Virtual file offset of the start of the chunk
  std::uint64_t chunk_end; //Virtual file offset of the end of the chunk
};

class Interval{
  std::uint64_t first_record_offset;
};

class Bin{
public:
  std::uint32_t bin_id;
  std::uint32_t v_offset_id;
  std::vector<Chunk> chunks;
};


class Index{
  std::vector<Bin> bins;
  std::vector<Interval> intervals_16kb;
  inline int reg2bins(int64_t beg, int64_t end, const std::int32_t min_interval_bits,const std::int32_t interval_idx_depth)
  {
    int l, t, n, s = min_interval_bits + interval_idx_depth*3;
    for (--end, l = n = t = 0; l <= interval_idx_depth; s -= 3, t += 1<<l*3, ++l) {
      int b = t + (beg>>s), e = t + (end>>s), i;
      for (i = b; i <= e; ++i) bins[n++].bin_id = i;
    }
    return n;
  }
  int reg2bin(int beg, int end)
  {
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
  }
  int reg2bins(int beg, int end, uint16_t lst[37449])
  {
    int i = 0, k;
    --end;
    lst[i++] = 0;
    for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) lst[i++] = k;
    for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) lst[i++] = k;
    for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) lst[i++] = k;
    for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) lst[i++] = k;
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) lst[i++] = k;
    return i;
  }
  inline int reg2bins(int64_t beg, int64_t end)
  {
    --end;
    constexpr std::int32_t min_interval_bits = 14; // bits for the minimal interval
    constexpr std::int32_t interval_idx_depth = 5; //Depth of the binning index
    int s = min_interval_bits + interval_idx_depth*3;
    int n=0;
    for (int l, t = 0; l <= interval_idx_depth; s -= 3) {
      int b = t + (beg>>s);
      int e = t + (end>>s), i;
      for (i = b; i <= e; ++i) bins[n++].bin_id = i;
      t += 1<<l*3;
      ++l;
    }
    return n;
  }
};

class CSIIndex{
  static constexpr std::int32_t min_interval_bits = 14; // bits for the minimal interval
  static constexpr std::int32_t interval_idx_depth = 5; //Depth of the binning index
  std::int32_t n_ref;
  std::int32_t col_seq;
  std::int32_t col_beg;
  std::int32_t col_end;
  std::int32_t meta;
  std::int32_t skip;
  std::int32_t l_nm;
  std::vector<char> seqname_buff;
  std::vector<Index> indices;
  /* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
  inline int reg2bin(int64_t beg, int64_t end)
  {
  int l, s = min_interval_bits;
  int t = ((1<<interval_idx_depth*3) - 1) / 7;
  for (--end, l = interval_idx_depth; l > 0; --l, s += 3, t -= 1<<l*3)
    if (beg>>s == end>>s)
      return t + (beg>>s);
  return 0;
}
  /* calculate the list of bins that may overlap with region [beg,end) (zero-based) */

/* calculate maximum bin number -- valid bin numbers range within [0,bin_limit) */
inline int bin_limit(int min_shift, int depth)
{
  return ((1 << (depth+1)*3) - 1) / 7;
}

};



class TabixIndex{
  std::int32_t n_ref;
  std::int32_t col_seq;
  std::int32_t col_beg;
  std::int32_t col_end;
  std::int32_t meta;
  std::int32_t skip;
  std::int32_t l_nm;
  std::vector<char> seqname_buff;
  std::vector<Index> indices;
};

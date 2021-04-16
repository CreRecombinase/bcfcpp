#pragma once

#include <cstdint>
#include <vector>
#include <fstream>


enum class endian
{
#ifdef _WIN32
    little = 0,
    big    = 1,
    native = little
#else
    little = __ORDER_LITTLE_ENDIAN__,
    big    = __ORDER_BIG_ENDIAN__,
    native = __BYTE_ORDER__
#endif
};

/*
  BGZF files support random access through the BAM file index. To achieve this,
  the BAM file index uses virtual file offsets into the BGZF file. Each virtual
  file offset is an unsigned 64-bit integer, defined as: coffset<<16|uoffset,
  where coffset is an unsigned byte offset into the BGZF file to the beginning
  of a BGZF block, and uoffset is an unsigned byte offset into the uncompressed
  data stream represented by that BGZF block. Virtual file offsets can be
  compared, but subtraction between virtual file offsets and addition between a
  virtual offset and an integer are both disallowed.
*/
struct __attribute__((__packed__)) voffset{
  std::uint64_t coffset : 48;
  std::uint64_t uoffset : 16;
};

class Chunk{
  std::uint64_t chunk_start; //Virtual file offset of the start of the chunk
  std::uint64_t chunk_end; //Virtual file offset of the end of the chunk
};

class Interval{
  std::uint64_t first_record_offset;
};

/*

The UCSC binning scheme was suggested by Richard Durbin and Lincoln Stein and is explained in Kent
et al.34 In this scheme, each bin represents a contiguous genomic region which is either fully contained in
or non-overlapping with another bin; each alignment is associated with a bin which represents the smallest
region containing the entire alignment. The binning scheme is essentially a representation of R-tree. A
distinct bin uniquely corresponds to a distinct internal node in a R-tree. Bin A is a child of Bin B if the
region represented by A is contained in B.
To find the alignments that overlap a specified region, we need to get the bins that overlap the region,
and then test each alignment in the bins to check overlap. To quickly find alignments associated with a
specified bin, we can keep in the index the start file offsets of chunks of alignments which all have the bin.
As alignments are sorted by the leftmost coordinates, alignments having the same bin tend to be clustered
together on the disk and therefore usually a bin is only associated with a few chunks. Traversing all the
alignments having the same bin usually needs a few seek calls. Given the set of bins that overlap the specified
region, we can visit alignments in the order of their leftmost coordinates and stop seeking the rest when an
alignment falls outside the required region. This strategy saves half of the seek calls in average.
In the BAI format, each bin may span 229, 226, 223, 220, 217 or 214 bp. Bin 0 spans a 512Mbp region, bins
1–8 span 64Mbp, 9–72 8Mbp, 73–584 1Mbp, 585–4680 128Kbp, and bins 4681–37448 span 16Kbp regions.
This implies that this index format does not support reference chromosome sequences longer than 229 − 1.
The CSI format generalises the sizes of the bins, and supports reference sequences of the same length as
are supported by SAM and BAM.


 */
class Bin{
public:
  std::uint32_t distinct_bin;
  std::uint32_t v_offset;
  std::vector<Chunk> chunks;
  Bin(std::ifstream &ifs){
    ifs.read((char*) &distinct_bin,sizeof(distinct_bin));
    ifs.read((char*) &v_offset,sizeof(v_offset));
    int32_t n_chunks;
    ifs.read((char*) &n_chunks,sizeof(n_chunks));
    chunks.resize(n_chunks);
    ifs.read((char*) chunks.data(),chunks.size()*sizeof(Chunk));
    //(char[sizeof(voffset)])"bla";
  }

};


class Index{
public:
  std::vector<Bin> bins;
  std::vector<Interval> intervals_16kb;
  Index(std::ifstream& ifs){
    std::int32_t n_bins;
    ifs.read((char*) &n_bins,sizeof(n_bins));
    if(n_bins<0)
      throw std::runtime_error("Index constructor encountered negative bin count:'"+std::to_string(n_bins));
    std::uint32_t distinct_bin;
    std::uint32_t distinct_bin;



  }
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

struct hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
  //    uint64_t n_no_coor;
    // bidx_t **bidx;
    // lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    int tbi_n, last_tbi_tid;

};


class idx_cache_t{
public:
  uint32_t last_bin, save_bin;
  std::int64_t last_coor;
  int last_tid, save_tid, finished;
  uint64_t last_off, save_off;
  uint64_t off_beg, off_end;
  uint64_t n_mapped, n_unmapped;
  idx_cache_t(uint64_t offset0):
    last_bin(0xffffffffu),
    save_bin(0xffffffffu),
    last_coor(0xffffffffu),
    last_tid(-1),
    save_tid(-1),
    finished(0),
    last_off(offset0),
    save_off(offset0),
    off_beg(offset0),
    off_end(offset0){};

}; // keep internal states


struct __attribute__((__packed__)) CSI_header{
  char magic[4];
  int32_t min_shift;
  int32_t depth;
  int32_t n_aux;
};

inline CSI_header read_CSI_header(std::ifstream& ifs){
  CSI_header csif;
  ifs.read((char*) &csif,sizeof(csif));
  if(memcmp(csif.magic, "CSI\1", 4) != 0)
    throw std::runtime_error("Incorrect magic str:'"+std::string(csif.magic)+std::string("CSI\1"));
  // if(x[0]!=min_interval_bits)      throw std::runtime_error("Only default min_interval_bits currently supported:'"+std::to_string(x[0])+" vs "+std::to_string(min_interval_bits));
  //   if(x[1]!=interval_idx_depth)      throw std::runtime_error("Only default interval_idx_depth currently supported:'"+std::to_string(x[1])+" vs "+std::to_string(interval_idx_depth));
  return csif;
}

class CSIIndex{
  static constexpr std::int32_t min_interval_bits = 14; // bits for the minimal interval
  static constexpr std::int32_t interval_idx_depth = 5; //Depth of the binning index
public:
  std::int32_t n_ref;
  std::int32_t col_seq;
  std::int32_t col_beg;
  std::int32_t col_end;
  std::int32_t meta;
  std::int32_t skip;
  std::int32_t l_nm;
  std::vector<char> aux_buff;
  //  std::vector<char> seqname_buff;
  std::vector<Index> indices;
  mutable idx_cache_t z;
  inline CSIIndex(std::ifstream & ifs):z(0){

    auto csif = read_CSI_header(ifs);
    aux_buff.resize(csif.n_aux);
    int32_t n_sequences;
    ifs.read((char*) &n_sequences,sizeof(n_sequences));
    if (n_sequences > INT32_MAX) throw std::overflow_error("overflow in number of possible sequences");
    indices.resize(n_sequences);

    z.save_tid = z.last_tid = -1;
    z.save_bin = z.last_bin = 0xffffffffu;
    z.save_off = z.last_off = z.off_beg = z.off_end = 0;
    z.last_coor = 0xffffffffu;
  }


    // if ((idx = hts_idx_init(n, HTS_FMT_CSI, 0, x[0], x[1])) == NULL) goto fail;
    // idx->l_meta = x[2];
    // idx->meta = meta;
    // meta = NULL;
    // if (idx_read_core(idx, fp, HTS_FMT_CSI) < 0) goto fail;

  inline constexpr std::int32_t n_bins(){
    return ((1<<(3 * interval_idx_depth + 3)) - 1) / 7;
  };
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

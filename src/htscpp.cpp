#include "htslib/khash_str2int.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/synced_bcf_reader.h"
#include <vector>


int read_line_bcf(std::vector<std::byte> &out_vec,BGZF* fp){
  uint8_t x[32];

  ssize_t ret;
  uint32_t shared_len, indiv_len;
  if ((ret = bgzf_read(fp, x, 32)) != 32) {
    if (ret == 0) return -1;
    return -2;
  }
  bcf_clear1(v);
  shared_len = le_to_u32(x);
  if (shared_len < 24) return -2;
  shared_len -= 24; // to exclude six 32-bit integers
  if (ks_resize(&v->shared, shared_len) != 0) return -2;
  indiv_len = le_to_u32(x + 4);
  if (ks_resize(&v->indiv, indiv_len) != 0) return -2;
  v->rid  = le_to_i32(x + 8);
  v->pos  = le_to_u32(x + 12);
  v->rlen = le_to_i32(x + 16);
  v->qual = le_to_float(x + 20);
  v->n_info = le_to_u16(x + 24);
  v->n_allele = le_to_u16(x + 26);
  v->n_sample = le_to_u32(x + 28) & 0xffffff;
  v->n_fmt = x[31];
  v->shared.l = shared_len;
  v->indiv.l = indiv_len;
  // silent fix of broken BCFs produced by earlier versions of bcf_subset, prior to and including bd6ed8b4
  if ( (!v->indiv.l || !v->n_sample) && v->n_fmt ) v->n_fmt = 0;

  if (bgzf_read(fp, v->shared.s, v->shared.l) != v->shared.l) return -2;
  if (bgzf_read(fp, v->indiv.s, v->indiv.l) != v->indiv.l) return -2;
  return 0;


}

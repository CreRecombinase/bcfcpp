#pragma once

#include "htslib/khash_str2int.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/synced_bcf_reader.h"
#include "range/v3/view/counted.hpp"
#include "range/v3/view/transform.hpp"
#include <math.h>

#include <iostream>
#include <istream>
#include <optional>
#include <string>
#include <variant>
#include <stdexcept>


#include <range/v3/range_fwd.hpp>

#include <range/v3/iterator/default_sentinel.hpp>
#include <range/v3/utility/static_const.hpp>
#include <range/v3/view/facade.hpp>
#include <range/v3/view/span.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/split_when.hpp>
#include <range/v3/view/for_each.hpp>


#include <range/v3/view/delimit.hpp>


#include <range/v3/detail/prologue.hpp>
#include <string_view>

template <class To, class From>
typename std::enable_if_t<
    sizeof(To) == sizeof(From) &&
    std::is_trivially_copyable_v<From> &&
    std::is_trivially_copyable_v<To>,
    To>
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
    static_assert(std::is_trivially_constructible_v<To>,
        "This implementation additionally requires destination type to be trivially constructible");

    To dst;
    std::memcpy(&dst, &src, sizeof(To));
    return dst;
}

class BCFLine{
public:
  bcf1_t *line;
  BCFLine():line(bcf_init1()){};
  ~BCFLine(){
    bcf_destroy1(line);
  };
};


template <typename T>
T convert_le(uint8_t *buf) {

  return reinterpret_cast<T>(buf);
}

template <> int8_t convert_le<int8_t>(uint8_t *buf) { return le_to_i8(buf); }

template <> int16_t convert_le<int16_t>(uint8_t *buf) { return le_to_i16(buf); }

template<>
int32_t convert_le<int32_t>(uint8_t *buf){
  return le_to_i32(buf);
}

template<>
float convert_le<float>(uint8_t *buf){
  return bit_cast<float>(le_to_u32(buf));
}

template <typename T> struct bcf_vector {
  static const T end;
  static const T missing;
};

template <> struct bcf_vector<int8_t> {
  static const int8_t end=bcf_int8_vector_end;
  static const int8_t missing=bcf_int8_missing;
};

template <> struct bcf_vector<int16_t> {
  static const int16_t end=bcf_int16_vector_end;
  static const int16_t missing=bcf_int16_missing;
};


template <> struct bcf_vector<int32_t> {
  static const int32_t end=bcf_int32_vector_end;
  static const int32_t missing=bcf_int32_missing;
};


// template <> struct bcf_vector<float> {
//   static const float  end=bcf_float_vector_end;
//   static const float missing=bcf_float_missing;
// };

//float bcf_float_missing    = memcpy(reinterpret_cast<float>(0x7F800001u);

// HTSLIB_EXPORT
// uint32_t bcf_float_vector_end = 0x7F800002;

//template <typename T>
inline bool is_vector_end(int8_t x) { return x == bcf_vector<int8_t>::end; }
inline bool is_vector_end(int16_t x){
  return x==bcf_vector<int16_t>::end;
}
inline bool is_vector_end(int32_t x) { return x == bcf_vector<int32_t>::end; }
inline bool is_vector_end(float x) {
  const auto t = std::numeric_limits<float>::signaling_NaN();
  return std::memcmp((void*)&x,(void*)&t,sizeof(float));
}

template<typename T>
class LineFMT{
public:
  bcf_fmt_t *fmt;
  //  ranges::span<T> data;
  LineFMT(bcf_fmt_t *fmt_):fmt(fmt_){
  };

  //returns a range of ranges with the outer-most range being a sample-level range
  auto sample_range(int n_samples)const {
    ranges::span<unsigned char*> data(&fmt->p,fmt->n);

    return ranges::views::chunk(data,sizeof(T)) |
      ranges::views::for_each([](unsigned char* x){
        return ranges::yield(convert_le<T>(x));
      }) | ranges::views::chunk(n_samples) | ranges::views::delimit([](auto &r){
        return is_vector_end(r);
      });
  };
};


using FMT_v = std::variant <std::monostate,
                            LineFMT<std::int8_t>,
                            LineFMT<std::int16_t>,
                            LineFMT<std::int32_t>,
                            LineFMT<std::int64_t>,
                            LineFMT<float>,
                            LineFMT<char> >;


class HTSFile{
public:
  htsFile* file;
  HTSFile(std::string_view file_name,std::string_view args):file(hts_open(file_name.data(),args.data())){};
  ~HTSFile(){
    hts_close(file);
  };
};



class BCF_ID_Tag {
public:
  std::string_view key;
  bcf_idinfo_t* val;
  BCF_ID_Tag(char* key_,bcf_idinfo_t *val_):key(key_),val(val_){
  };


  // int id_int(const BCFHeader& hdr,) const {
  //   bcf_hdr_id2int(
};



template<int HLType>
class HeaderLine{
public:
  const bcf_idpair_t* id;
  HeaderLine(const bcf_idpair_t* id_):id(id_){};
  std::string_view name() const {return id->key;}
  int type_length_i() const{
    if constexpr(HLType==0)
      return 0;
    auto idx= static_cast<int>(id->val->info[HLType]>>8 & 0xf);
    //constexpr std::string_view sv[] {"Fixed","Variable","Allele","Genotype","Region"};
    //    return sv[idx];

    return idx;
  }
  std::string_view type_length_s() const{
    constexpr std::string_view sv[] {"Fixed","Variable","Allele","Genotype","Region"};
    return sv[type_length_i()];
  }
  int type_number() const{
    return static_cast<int>(id->val->info[HLType]>>12);
  }
  int type_type_i() const{
    if constexpr(HLType==0)
      return 0;
    auto idx = static_cast<int>(id->val->info[HLType]>>4 & 0xf);
    return idx;
  }
  std::string_view type_type_s() const{
    constexpr std::string_view sv[] {"BOOLEAN","INT","REAL","VARCHAR"};
    return sv[type_type_i()];
  }
  std::string db_create_name() const{
    return name()+" "+type_type_s();
  }
  constexpr int HL_type() const{
    return HLType;
  }
  constexpr std::string_view HL_type_s() const{
    constexpr std::string_view sv[] {"Filter","Info","Format"};
    return sv[HL_type()];
  }

};

using vl_var= std::variant<HeaderLine<0>,HeaderLine<1>,HeaderLine<2> >;

inline vl_var idpair2HL(const bcf_idpair_t *idp){
  if(static_cast<uint32_t>(idp->val->info[0] & 0xf)==0){
    return HeaderLine<0>(idp);
  }
  if(static_cast<uint32_t>(idp->val->info[1] & 0xf)==1){
    return HeaderLine<1>(idp);
  }
  if(static_cast<uint32_t>(idp->val->info[2] & 0xf)!=2){
    throw std::runtime_error("bcf_idpair_t is of invalid type!");
  }
  return HeaderLine<2>(idp);
}

template <typename T>
class BCFInfo{
public:
  bcf_info_t* info;
};




class BCFHeader{
public:
  bcf_hdr_t * header;
  mutable std::optional<int> gt_id;
  BCFHeader(HTSFile &file):header(bcf_hdr_read(file.file)){}
  ~BCFHeader(){
    bcf_hdr_destroy(header);
  }
  std::optional<int> get_GT_id() const{
    if(!gt_id){
      int tag_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");
      if (bcf_hdr_idinfo_exists(header,BCF_HL_FMT,tag_id) )
        gt_id=tag_id;    // no such FORMAT field in the header
    }
    return gt_id;
  }
  // template<int which>
  // int get_idint(const char* id); const{
  //   static_assert(which>0 ,"get_idint must be with nonnegative value for 'which'");
  //   static_assert(which<=2 ,"get_idint<which> must have 'which' <=2");
  //   khint_t k;
  //   vdict_t *d = (vdict_t*)h->dict[which];
  //   k = kh_get(vdict, d, id);
  //   return k == kh_end(d)? -1 : kh_val(d, k).id;
  // }


  int32_t num_samples()const{
    return header->n[2];
  }
  int32_t num_contigs()const {
    return header->n[1];
  }
  int32_t num_ID()const {
    return header->n[0];
  }
  int id2length(int idx) const{
    return bcf_hdr_id2length(header,BCF_HL_INFO,idx);
  }
  auto get_IDs() const {
    const ::ranges::span<bcf_idpair_t> counted_r(header->id[0],header->n[0]);
    auto transformed_r=ranges::views::transform(counted_r,[](const bcf_idpair_t & idp){
      return idpair2HL(&idp);
    });
    return transformed_r;
  }
  auto get_INFOs() const {
    const ::ranges::span<bcf_idpair_t> counted_r(header->id[0],header->n[0]);
    auto transformed_r=ranges::views::transform(counted_r,[](const bcf_idpair_t & idp){
      return idpair2HL(&idp);
    });
    auto filtered_r = ranges::views::filter(transformed_r,[](const vl_var &v) {
      return std::visit([](auto && vlvv) -> bool { return vlvv.HL_type()==1;},v);
    });
    return ranges::views::transform(filtered_r,[](const vl_var &v){
      return std::get<HeaderLine<1>>(v);
    });
  }
  ranges::subrange<bcf_hrec_t **, bcf_hrec_t **, ranges::subrange_kind::sized>
  get_hrecs() const{
    return ranges::views::counted(header->hrec,header->nhrec);
  }
  int32_t num_header_records() const{
    return header->nhrec;
  }
};


class BCFFile{
public:
  HTSFile file;
  BCFHeader header;

  BCFFile(std::string_view file_name, std::string_view args="r"):file(file_name,args),header(file){};

  
};


template<int Fields=BCF_UN_ALL>
class UnpackedBCFLine:public BCFLine{
public:
  UnpackedBCFLine(){
    //    bcf_unpack(line,Fields);
  };
  std::optional<std::string_view> get_ID(){
    if(line->d.id[0]=='.')
      return std::nullopt;
    return std::string_view(line->d.id);
  }
  int get_line_id(const int header_id) const{
    const bcf_fmt_t *fmtb = line->d.fmt;
    const bcf_fmt_t* fmte = line->d.fmt+line->n_fmt;
    auto fr = std::find(fmtb,fmte,[&header_id](const bcf_fmt_t &fmt){
      return fmt.id==header_id;
    });
    return fr-fmtb;
  }
  FMT_v get_FMT(int i){
    auto fmt = &line->d.fmt[i];
    if(!fmt->p){
      return std::monostate{};
    }
    switch (fmt->type) {
    case BCF_BT_INT8:
      return LineFMT<std::int8_t>(fmt);
    case BCF_BT_INT16:
      return LineFMT<std::int16_t>(fmt);
    case BCF_BT_INT32:
      return LineFMT<std::int32_t>(fmt);
    case BCF_BT_INT64:
      return LineFMT<std::int64_t>(fmt);
    case BCF_BT_FLOAT:
      return LineFMT<float>(fmt);
    default:
      throw std::invalid_argument("Unexpected type "+std::to_string(fmt->type));
    }
  }


  int get_GT(const BCFFile & file, std::string &dest) const{
    auto tag_id = file.header.get_GT_id();
    if(!tag_id.has_value()){
      return -1;
    }

    int i= get_line_id(*tag_id);                    // the tag is not present in this record
    bcf_fmt_t *fmt = &line->d.fmt[i];
    if ( !fmt->p ) return -3;                                      // the tag was marked for removal


    int n = fmt->n*file.header.num_samples();
    char* new_p = bit_cast<char*>(&(fmt->p));
    dest=std::string_view(new_p,n);
    return n;

  }

  template<typename T>
  int copy_GT(const BCFFile & file, std::vector<T> &dest) const{
    auto tag_id = file.header.get_GT_id();
    if(!tag_id.has_value()){
      return -1;
    }

    int i= get_line_id(*tag_id);
    auto fmtv = get_FMT(i);



    // Make sure the buffer is big enough
    int nsmpl = file.header.num_samples();
    dest.reserve(nsmpl);


    #define BRANCH(type_t, convert, is_missing, is_vector_end, set_missing, set_vector_end, set_regular, out_type_t) { \
        out_type_t *tmp = (out_type_t *) *dst; \
        uint8_t *fmt_p = fmt->p; \
        for (i=0; i<nsmpl; i++) \
        { \
            for (j=0; j<fmt->n; j++) \
            { \
                type_t p = convert(fmt_p + j * sizeof(type_t)); \
                if ( is_missing ) set_missing; \
                else if ( is_vector_end ) { set_vector_end; break; } \
                else set_regular; \
                tmp++; \
            } \
            for (; j<fmt->n; j++) { set_vector_end; tmp++; } \
            fmt_p += fmt->size; \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8, p==bcf_int8_missing,  p==bcf_int8_vector_end,  *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, *tmp=p, int32_t); break;
        case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, p==bcf_int16_missing, p==bcf_int16_vector_end, *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, *tmp=p, int32_t); break;
        case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, p==bcf_int32_missing, p==bcf_int32_vector_end, *tmp=bcf_int32_missing, *tmp=bcf_int32_vector_end, *tmp=p, int32_t); break;
        case BCF_BT_FLOAT: BRANCH(uint32_t, le_to_u32, p==bcf_float_missing, p==bcf_float_vector_end, bcf_float_set_missing(*tmp), bcf_float_set_vector_end(*tmp), bcf_float_set(tmp, p), float); break;
        default: hts_log_error("Unexpected type %d at %s:%" PRIhts_pos, fmt->type, bcf_seqname_safe(hdr,line), line->pos+1); exit(1);
    }
    #undef BRANCH
    return nsmpl*fmt->n;
}



}



};






inline int getline_bcf(BCFFile &bcf_file,BCFLine &line){
  return bcf_read((bcf_file.file.file), (bcf_file.header.header), (line.line));
}

template <int Fields=BCF_UN_ALL>
inline int getline_bcf(BCFFile &bcf_file,UnpackedBCFLine<Fields> &line){
  int ret =bcf_read((bcf_file.file.file), (bcf_file.header.header), (line.line));
  if(ret==0){
    bcf_unpack(line.line,Fields);
  }
  return ret;
}






namespace ranges
{
    /// \addtogroup group-views
    /// @{
  template<int Fields=BCF_UN_ALL>
    struct getlines_hts_view : view_facade<getlines_hts_view<Fields>, unknown>
    {
    private:
      friend range_access;
      BCFFile * sin_;
      UnpackedBCFLine<Fields> str_;
      struct cursor
      {
      private:
        friend range_access;
        using single_pass = std::true_type;
        getlines_hts_view * rng_ = nullptr;

      public:
        cursor() = default;
        explicit cursor(getlines_hts_view * rng)
          : rng_(rng)
        {}
        void next()
        {
          rng_->next();
        }
        UnpackedBCFLine<Fields> & read() const noexcept
        {
          return rng_->str_;
        }
        bool equal(default_sentinel_t) const
        {
          return !rng_->sin_;
        }
        bool equal(cursor that) const
        {
          return !rng_->sin_ == !that.rng_->sin_;
        }
      };
      void next()
      {
        if(getline_bcf(*sin_,str_)==-1){
          sin_ = nullptr;
        }else{

        }
      }
      cursor begin_cursor()
      {
        return cursor{this};
      }

    public:
      getlines_hts_view() = default;
      getlines_hts_view(BCFFile & sin)
        : sin_(&sin)
        , str_{}
      {
        this->next(); // prime the pump
      }
      BCFLine & cached() noexcept
      {
        return str_;
      }
    };

template <int Fields=BCF_UN_ALL>
  struct getlines_hts_fn
  {
    getlines_hts_view<Fields> operator()(BCFFile & sin, int which=BCF_UN_STR) const
    {
      return getlines_hts_view<Fields>{sin, which};
    }
  };

  RANGES_INLINE_VARIABLE(getlines_hts_fn, getlines_hts)
  /// @}
} // namespace ranges

#include <range/v3/detail/epilogue.hpp>

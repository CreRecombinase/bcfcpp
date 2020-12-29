#pragma once


#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/synced_bcf_reader.h"
#include "range/v3/view/counted.hpp"
#include "range/v3/view/transform.hpp"
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

#include <range/v3/detail/prologue.hpp>
#include <string_view>

class BCFLine{
public:
  bcf1_t *line;
  BCFLine():line(bcf_init1()){};
  ~BCFLine(){
    bcf_destroy1(line);
  };
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
};


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
  BCFHeader(HTSFile &file):header(bcf_hdr_read(file.file)){}
  ~BCFHeader(){
    bcf_hdr_destroy(header);
  }
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

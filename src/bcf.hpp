#pragma once

#include <string>
#include <optional>
#include <vector>
#include <array>
#include <fstream>

struct __attribute__((__packed__)) BCF_Header{
  char magic[3];
  char major_version;
  char minor_version;
  std::uint32_t header_size;
};

inline bool is_valid_header(const BCF_Header & h){
  if((h.magic[0]!='B') or (h.magic[1]!='C') or (h.magic[2]!='F') or h.major_version!=2 )
    return false;
  return true;
}

struct __attribute__((__packed__)) BCFBuff{
  std::uint32_t l_shared;
  std::uint32_t l_indiv;
  int32_t chrom;
  int32_t pos;
  int32_t rlen;
  float qual;
  int32_t n_info: 16;
  int32_t n_allele: 16;
  uint32_t n_sample : 24;
  uint32_t n_fmt : 8;
};



class BCF{
public:
  const char magic_str[5] {'B','C','F','\2','\2'};
  inline  const char* magic_string(){
    return magic_str;
  };
  std::string header_text;
  BCF(std::ifstream & ifs){
    std::array<char,5> check_magic;
    BCF_Header hdr;
    ifs.read((char*) &hdr, sizeof(hdr));
    if(!is_valid_header(hdr))
      throw std::runtime_error("Incorrect magic str:'"+std::string(hdr.magic)+std::to_string(hdr.major_version)+"."+std::to_string(hdr.minor_version)+"'vs'"+std::string(magic_str)+"'");

    header_text.resize(hdr.header_size);
    ifs.read(header_text.data(),header_text.size());
  }

};


class BCFRec{

  int32_t chrom;
  int32_t pos;
  std::optional<float> qual;
  std::optional<std::string> id;
  std::vector<int> filtervec;
};

class BCFReader {
  BCFRec bcfb;
  std::ifstream& ifs;
  std::vector<std::byte> shared_buff;
  std::vector<std::byte> indiv_buff;
  int getline_bcf_record(std::ifstream & ifs,BCFBuff &bcfb){
    if(!ifs)
      return 0;
    ifs.read((char*) &bcfb,sizeof(bcfb));

    shared_buff.resize(bcfb.l_shared);
    indiv_buff.resize(bcfb.l_indiv);
    if(!ifs.read((char*)shared_buff.data(),shared_buff.size()))
      return -1;
    ifs.read((char*)indiv_buff.data(),indiv_buff.size());
    return 1;
  }










};

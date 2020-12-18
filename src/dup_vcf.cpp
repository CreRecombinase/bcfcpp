#include "cxxopts.hpp"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/synced_bcf_reader.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <iostream>
#include <range/v3/view/all.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/indices.hpp>
#include <string_view>
#include <optional>



//#include <range/v3/all.hpp>

void find_and_increment(std::string_view sv,std::unordered_map<size_t,char>& ssv,int j){
  auto svf = std::hash<std::string_view>{}(sv);

  auto fsf = ssv.find(svf);

  if(fsf==ssv.end()){
    auto [new_it,check]=ssv.insert({svf,0});
    fsf=new_it;
  }else{
    if(fsf->second==0){
      std::cout<<sv<<std::endl;
      fsf->second=1;
    }
  } 
}



int main(int argc, char** argv){

  cxxopts::Options options("dup_vcf", "Check a VCF/BCF for duplicate IDS");

  options.add_options()
    ("f,file", "VCF/BCF input file", cxxopts::value<std::vector<std::string>>())
    ("s,sync", "whether to use the synced reader", cxxopts::value<bool>()->default_value("false"))
    ("f,file", "VCF/BCF input file", cxxopts::value<std::vector<std::string>>())
    ("h,help", "Print usage");
  auto result = options.parse(argc, argv);
   if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
   std::vector<std::string> vcf_files(result.count("file"));
   if(vcf_files.size()==0){
      std::cerr<<"must specify at least 1 bcf file"<<std::endl;
      return 1;
    }
   vcf_files=result["file"].as<std::vector<std::string>>();
   if(result.count("sync")>0){
     std::unique_ptr<bcf_srs_t,decltype(&bcf_sr_destroy)> sr(bcf_sr_init(),&bcf_sr_destroy);
     std::for_each(vcf_files.begin(),vcf_files.end(),[&sr](const auto &i){
       bcf_sr_add_reader(sr.get(),i.c_str());
     });
     std::unordered_map<size_t,char> all_ids;
     auto enumerated_vcfs = ranges::views::enumerate(vcf_files);
     std::string_view lds;
     int j=0;
     while ( bcf_sr_next_line(sr.get()) ){
       for(auto it = enumerated_vcfs.begin(); it != enumerated_vcfs.end(); ++it)
         {
           const auto idx = std::get<0>(*it);
           bcf1_t *line = bcf_sr_get_line(sr,idx);
           int bcfr=bcf_unpack(line,BCF_UN_STR);
           if(line->d.id[0]!='.'){
             lds = line->d.id;
             find_and_increment(lds,all_ids,j++);
           }
         }
     }
   }else{
     auto hts_rng = ranges::views::all(vcf_files) | ranges::views::transform([](const std::string& str){
       return std::unique_ptr<htsFile,decltype(&hts_close)> (hts_open(str.c_str(),'r'),&hts_close);
     });

  return 0;
}

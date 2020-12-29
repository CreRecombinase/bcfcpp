#include "cxxopts.hpp"
#include "bcfcpp.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <iostream>
#include <range/v3/view/all.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/indices.hpp>
#include <range/v3/core.hpp>
#include <range/v3/algorithm/for_each.hpp>
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
   std::unordered_map<size_t,char> all_ids;
   if(result.count("sync")>0){
     std::unique_ptr<bcf_srs_t,decltype(&bcf_sr_destroy)> sr(bcf_sr_init(),&bcf_sr_destroy);
     std::for_each(vcf_files.begin(),vcf_files.end(),[&sr](const auto &i){
       bcf_sr_add_reader(sr.get(),i.c_str());
     });

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
     //     std::unique_ptr<bcf_srs_t,decltype(&bcf_sr_destroy)> sr(bcf_sr_init(),&bcf_sr_destroy);
     int j=0;
     for(auto it = vcf_files.begin(); it !=vcf_files.end(); it++){
       BCFFile bcf(std::string_view{*it},std::string_view{"r"});
       auto mrng = ranges::getlines_hts_view<BCF_UN_STR>(bcf);
       //for( auto line = mrng.begin(); line!=mrng.end(); line++){
       ranges::for_each(mrng,[&all_ids,&j](UnpackedBCFLine<BCF_UN_STR>& line) mutable{
         if(auto li =line.get_ID())
           find_and_increment(*li,all_ids,j++);

       });

     }
   }
  return 0;
}

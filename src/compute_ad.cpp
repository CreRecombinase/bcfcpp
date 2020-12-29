#include "cxxopts.hpp"
#include "bcfcpp.hpp"
#include "range/v3/view/cache1.hpp"
#include "range/v3/view/join.hpp"

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

  cxxopts::Options options("header_info", "Get the schema of a BCF/VCF header");

  options.add_options()
    ("n,samplenum", "number of samples in the dataset", cxxopts::value<int>())
    ("h,help", "Print usage");
  auto result = options.parse(argc, argv);
  if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

  if(result.count("samplenum")==0){
    std::cerr<<"must specify sample size"<<std::endl;
    return 1;
  }
  int samplesize = result["samplenum"].as<int>();
  std::vector<double> means(samplesize);
  auto rli = ranges::views::getlines(std::cin) | ranges::views::split('\t') |ranges::views::transform([](auto &ir){
    ranges::views::enumerate(ir);
  });







     //     std::unique_ptr<bcf_srs_t,decltype(&bcf_sr_destroy)> sr(bcf_sr_init(),&bcf_sr_destroy);
  int j=0;
  for(auto it = vcf_files.begin(); it !=vcf_files.end(); it++){
    BCFFile bcf(std::string_view{*it},std::string_view{"r"});
    std::cout<<"Number of IDs: "<<bcf.header.num_ID()<<", Contigs: "<<bcf.header.num_contigs()<<", Samples: "<<bcf.header.num_samples()<<", num header records: "<<bcf.header.num_header_records()<<std::endl;

    auto    id_range = bcf.header.get_INFOs();

    auto inp_sr = ranges::views::transform(id_range,[](auto idp) ->std::string{
      return std::string(idp.name())+" "+std::string(idp.type_type_s());
    });
    std::string db_create = ranges::views::join(ranges::views::cache1(inp_sr),',') | ranges::to<std::string>();
    std::cout<<"CREATE TABLE info("<<db_create<<")"<<std::endl;




         // std::cout<<static_cast<uint32_t>(val->info[0] & 0xf)<<","<<static_cast<uint32_t>(val->info[1] & 0xf)<<","<<static_cast<uint32_t>(val->info[2] & 0xf)<<std::endl;

         // //         if(idp->type)
         // auto keys = ranges::views::counted(idp->keys,idp->nkeys);
         // auto vals = ranges::views::counted(idp->vals,idp->nkeys);
         // auto zkv = ranges::views::zip(keys,vals);
         // ranges::for_each(zkv,[](auto zkp){
         //   auto [keyv,valuev] = zkp;
         //   std::cout<<keyv<<" : "<<valuev<<std::endl;
         // //         std::cout<<bcf_hdr_id2length(bcf.header.header,
         // });
       //       });

     }

     return 0;
}

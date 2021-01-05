#include "cxxopts.hpp"
#include "bcfcpp.hpp"
#include "range/v3/view/cache1.hpp"
#include "range/v3/view/join.hpp"
#include "range/v3/view/zip.hpp"
#include "hts_bcf_getlines.hpp"
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


int main(int argc, char** argv){

  cxxopts::Options options("header_info", "Get the schema of a BCF/VCF header");

  options.add_options()
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
  //  std::vector<double> means(samplesize);
  // auto rli = ranges::getlines(std::cin) | ranges::views::split('\t') |ranges::views::transform([](auto &ir){
  //   ranges::views::enumerate(ir);
  // });







     //     std::unique_ptr<bcf_srs_t,decltype(&bcf_sr_destroy)> sr(bcf_sr_init(),&bcf_sr_destroy);
  int j=0;
  for(auto it = vcf_files.begin(); it !=vcf_files.end(); it++){
    BCFFile bcf(std::string_view{*it},std::string_view{"r"});
    const auto nsamples = bcf.header.num_samples();
    const auto tag_id = bcf.header.get_GT_id();
    if(!tag_id.has_value()){
      std::cerr<<"File must have gt field";
      return 1;
    }
    auto sample_v = bcf.header.view_samples();

    auto mrng = ranges::getlines_hts_view<BCF_UN_FMT>(bcf) | ranges::views::transform([=](UnpackedBCFLine<BCF_UN_FMT> &line) {
      return line.get_v_FMT(*tag_id);
    });
    auto i_mrng=ranges::views::enumerate(mrng);
    ranges::for_each(i_mrng,[&](auto && fmtvi){
      auto [i,fmtv] = fmtvi;
      std::cout<<"Variant "<<i<<std::endl;
      std::visit([&](auto &&arg){
        auto sr = arg.sample_range(nsamples);

        auto sample_v_r = ranges::views::enumerate(sr);
        int s=0;
        ranges::for_each(sample_v_r,[](const auto &ix) mutable{;
            auto [vi,x] = ix;
            std::cout<<"Sample: "<<vi<<":"<<std::endl;
            ranges::for_each(x,[](const auto &v){
              std::cout<<v<<",";
            });
            std::cout<<std::endl;
        });
      },fmtv);
    });
  }

  return 0;
}

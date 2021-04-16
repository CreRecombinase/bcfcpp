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



int main(int argc, char** argv){

  cxxopts::Options options("semi_join", "Subset A VCF/BCF file using one (or more) second VCF files");

  options.add_options()
    ("x,target", "VCF/BCF input to filter", cxxopts::value<std::string>())
    ("y,query", "VCF/BCF(s) to use for subsetting", cxxopts::value<std::vector<std::string>>())
    ("h,help", "Print usage");
  auto result = options.parse(argc, argv);
  if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

  if(result.count("target")==0){
    std::cerr<<"must specify target b/vcf file"<<std::endl;
    return 1;
  }
  std::string target_vcf_file = result["target"].as<std::string>();
  if(result.count("query")==0){
    std::cerr<<"must specify query b/vcf file"<<std::endl;
    return 1;
  }
  std::string query_vcf_file = result["target"].as<std::string>();
  //  std::vector<double> means(samplesize);
  // auto rli = ranges::getlines(std::cin) | ranges::views::split('\t') |ranges::views::transform([](auto &ir){
  //   ranges::views::enumerate(ir);
  // });

  BCFFile target_bcf(std::string_view{target_vcf_file},std::string_view{"r"});
  BCFFile query_bcf(std::string_view{query_vcf_file},std::string_view{"r"});

  BCFFile out_bcf(std::string_view{query_vcf_file},std::string_view{"w"});


}

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

  cxxopts::Options options("left_join", "Subset A VCF/BCF file using a second VCF file");

  options.add_options()
    ("a,target", "VCF/BCF input to filter", cxxopts::value<std::vector<std::string>>())
    ("b,query", "VCF/BCF sites at ", cxxopts::value<std::vector<std::string>>())
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


}

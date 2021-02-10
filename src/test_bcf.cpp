#include "bcf.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include "cxxopts.hpp"

int main(int argc, char** argv){

    cxxopts::Options options("test_bcf", "Read bcf without htslib");
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
  std::ifstream is(*(vcf_files.begin()),std::ios::binary);
  try
{

  BCF bcf(is);
  BCFBuff bcfb;
  is.read((char*) &bcfb,sizeof(bcfb));
  std::cout<<"l_shared:"<<bcfb.l_shared<<
    "\nl_indiv:"<<bcfb.l_indiv<<
    "\nchrom:"<<bcfb.chrom<<
    "\npos:"<<bcfb.pos<<
    "\nrlen:"<<bcfb.rlen<<
    "\nqual:"<<bcfb.qual<<
    "\nn_allele:"<<bcfb.n_allele<<
    "\nn_info:"<<bcfb.n_info<<
    "\nn_fmt:"<<bcfb.n_fmt<<
    "\nn_sample:"<<bcfb.n_sample<<std::endl;

 } catch (const std::exception& e) {
    std::cout << "main() failed to create C with: " << e.what();
  }
  return 0;






}

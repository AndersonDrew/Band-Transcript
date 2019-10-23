##Supporting Information 3

#include <Rcpp.h>
#include <cstdlib>
#include <string>
#include <vector>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
CharacterVector blastAnnotate(DataFrame blast, DataFrame gff) {
  Rcpp::IntegerVector blast_start = blast["sstart"];
  Rcpp::CharacterVector blast_chrom  = blast["sallseqid"];
  Rcpp::IntegerVector blast_end = blast["send"];
  
  Rcpp::IntegerVector gff_start = gff["start"];
  Rcpp::IntegerVector gff_end = gff["end"];
  Rcpp::CharacterVector gff_chrom = gff["seqname"];
  Rcpp::CharacterVector gff_attributes = gff["attribute"];
  
  std::vector<std::string> att_list;
  
  for(int i=0; i<blast_chrom.size();i++){
    att_list.push_back("");
    std::vector<std::string> attributes;
    for(int ii = 0; ii < gff_chrom.size(); ii++){
      if(blast_chrom[i] == gff_chrom[ii]){
        if(blast_start[i] >= gff_start[ii] && blast_end[i] <= gff_end[ii]){//then we keep the attributes
          attributes.push_back(Rcpp::as<std::string>(gff_attributes[ii]));
        }
      }
    }//done matching gff to blast
    if(attributes.size() > 1){
      for(std::size_t ii = 0; ii < attributes.size(); ii++){
        att_list[i] = att_list[i] + attributes[ii] + ";";
      }
    }else{
      att_list[i] = attributes[0];
    }
  }// done with for loop through blast_chroms
  CharacterVector IDs = wrap(att_list);
  return IDs;
}




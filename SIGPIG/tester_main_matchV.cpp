#include <iostream>
#include <vector>
#include <cstring>
#include <assert.h>

#include "MatrixFunctions.hpp"
#include "Alignment.hpp"
#include "Alignment_V.hpp"
#include "AlignmentSettings.hpp"

//this is a testing module for testing match v functions

using namespace std;
int main()
{
  
  //testing align with constraints fixed left remove right errors
  cout<<"%%%%%%%%testing alignment function 1 fixed left remove"<<endl;
  double error_cost=4;
  /*  string seq1="CAAGCCCCTTCATCAGGAACTGGTCCGTCGGGTCCCGGCGACACGGGGGTCTCCACGAGAACCTCCTCCCACGATCCCCCTTCTGGCTACCCGGGAACCACCTCCGACTCCTCTGCCACTGGCGCCAGGGAACCGGGGTCTGCAGGTATCGCATCATCGGTTTGGGACCCCATTGGGGAGGTTTTTAGCACTAGGGCCGACTCACCAGAGAGCGTGTCATTATGTGACGGCACAGGAGCCGAGCGGGGTACATGAGACGCAACTATGGTGACTAACGGGATATCAATGAGCAGACAGAGA";
  string seq2=    "GACTTCTCTGCCACTGGTAACAGGGAACCGGGGTCTGTAGTTTTCGTAGT";
  */
  /*  string seq1="GACTCCTCTGCCACTGGCGCCAGGGAACCGGGGTCTGCAGGTATCGCATCATC";
  string seq2="GACTCCTCTGTCACTGGTCCCACGGTGCCGGGGTCTCTAGCTTCATGGTCATC";
  */
  string seq1("AGAGACAGACGAGTAACTATAGGGCAATCAGTGGTATCAACGCAGAGTACATGGGGCGAGCCGAGGACACGGCAGTGTATTACTGTGCGAGAGACCACTCAGCCGGGATCACGATTTTTGGAGGGGTTACCCCAGGGTTTGGCTACTACGCTATGGACGTCTGGGGCCAAGGGACCGCGGTCACCGTCTCCTCAGCCTCCACCAAGGGCCCATCGGTCTTCCCCCTAGCACCCTCCTCCAAGAGCACCTCTGGGGGCACAGCGGCCCTGGGCTGCCTGGTCAAGGACTACTTCCCCGAAC");
  //string seq2("CAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGTGCTGGATCCGCCAGCCCCTAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACAACAACCCGTCCCTCAAGAGTCGAGCCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCTGTGTATTACTGTGCGAGAGG");
  //string seq2("GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTAGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCGTCAGTAGCAATGAGATGAGCTGGATCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTGGTGGTAGCACATACTACGCAGACTCCAGGAAGGGCAGATTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTTCAAATGAACAACCTGAGAGCTGAGGGCACGGCCGCGTATTACTGTGCCAGATATA");
  //string seq2("GAAGTGCAGCTGGTGGAGTCTGGGGGAGTCGTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTGGGTCCGTCAAGCTCCGGGGAAGGGTCTGGAGTGGGTCTCTCTTATTAGTTGGGATGGTGGTAGCACCTACTATGCAGACTCTGTGAAGGGTCGATTCACCATCTCCAGAGACAACAGCAAAAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACCGCCTTGTATTACTGTGCAAAAGATA");
  string seq2("GAAGTGCAGCTGGTGGAGTCTGGGGGAGTCGTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTGGGTCCGTCAAGCTCCGGGGAAGGGTCTGGAGTGGGTCTCTCTTATTAGTTGGGATGGTGGTAGCACCTACTATGCAGACTCTGTGAAGGGTCGATTCACCATCTCCAGAGACAACAGCAAAAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACCGCCTTGTATTACTGTGCAAAAGATA");
  unsigned maximum_errors=AlignmentSettings::V_allowed_errors;
  unsigned n_errors=0;
  unsigned* error_positions=new unsigned [maximum_errors];
  unsigned alength=20;
  unsigned* align_positions=new unsigned[2];
  unsigned align_position_start;

  cout<<"#####################start doing the calculate score"<<endl;
  double score= CalculateScore
    (alength, NULL, 5,4);
  assert(score==alength-5*4);
  cout<<"score is correct"<<endl;
  
  cout<<"#####################start doing align remove both errors:"<<endl;
  alength= align_with_constraints_fixed_left_remove_both_errors
    ( seq1, seq2,  AlignmentSettings::V_allowed_errors ,error_cost,
     /*output*/ n_errors, error_positions,
      align_position_start);
  cout<<"\talength:"<<alength<<endl;
  cout<<"\tn_errors:"<<n_errors<<endl;
  cout<<"\talign_position_start:"<<align_position_start<<endl;
  cout<<"\terror_positions:";
  for(unsigned int i=0;i<n_errors;i++)
    {
      cout<<error_positions[i]<<",";
    }
  cout<<endl;

  cout<<"#####################start doing align fast no fix:"<<endl;
  alength=align_with_constraints_fast_no_fix
     (seq1, seq2, AlignmentSettings::V_allowed_errors ,
      AlignmentSettings::V_minimum_alignment_length, error_cost,
      /*output*/align_positions, n_errors, error_positions);
  cout<<"\talength:"<<alength<<endl;
  cout<<"\tn_errors:"<<n_errors<<endl;
  cout<<"\talign_positions:"<<align_positions[0]
      <<","<<align_positions[1]<<endl;
  cout<<"\terror_positions:";
  for(unsigned int i=0;i<n_errors;i++)
    {
      cout<<error_positions[i]<<",";
    }
  cout<<endl;
  cout<<"\tmin_deletion:"
      <<seq2.size()-(align_positions[1]+alength)<<endl;
  
  //testing match V
  

  //clear the memory
  delete[] error_positions;
  delete[] align_positions;
  cout<<"!!!!!DONE!!!!!!"<<endl;
  
  return 0;
}

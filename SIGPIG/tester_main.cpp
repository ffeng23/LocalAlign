#include <iostream>
#include <vector>
#include <cstring>

#include "MatrixFunctions.hpp"
#include "Alignment.hpp"
//this is a testing module for testing various functions

using namespace std;
int main()
{
  //first testing matrixFunctions
  cout<<"Testing max and min functions:"<<endl;
  vector<unsigned> vec_int;
  vec_int.push_back(1);
  vec_int.push_back(2);
  vec_int.push_back(3);
  vec_int.push_back(0);
  
  vector<double> vec_double;
  vec_double.push_back(1);
  vec_double.push_back(2);
  vec_double.push_back(3);
  vec_double.push_back(0);
  
  unsigned array_int[]={1,2,3,0};
  double array_double[]={1,2,3,0};
  
  double array_double2[]={ 7, -13, 1, 3, 10, 5, 2, 4 };

  cout<<"return max : "<<max_mf(array_int, 4)<<endl;
  cout<<"return max : "<<max_mf(array_double, 4)<<endl;
  
  cout<<"return max : "<<max_mf(vec_int)<<endl;
  cout<<"return max : "<<max_mf(vec_double)<<endl;

  cout<<"return min : "<<min_mf(array_int, 4)<<endl;
  cout<<"return min : "<<min_mf(array_double, 4)<<endl;
  //const vector<unsigned> constVec=vec_int;
  cout<<"return min : "<<min_mf(vec_int)<<endl;
  cout<<"return min : "<<min_mf(vec_double)<<endl;

  cout<<"***********testing get median index***********"<<endl;
  cout<<"\tmedian:first:"<<0<<";last:"<<7<<";middle:"<<3<<";median:"<<GetMedianIndex(array_double2,0,7,3, NULL)<<endl;
  cout<<"\t\tpivot value:"<<array_double2[GetMedianIndex(array_double2,0,7,3, NULL)]<<endl;

  cout<<"*******Testing quick sort*******"<<endl;
  unsigned index[]={0,1,2,3,4,5,6,7};
  Print(array_double2, 8);
  QuickSort(array_double2, 0,7, index);
  Print(array_double2, 8);

  Print(index, 8);

  cout<<"****New array********************>>>>>>>>>>>>>"<<endl;
  double array_double3[]={ 1.0, 23.0, 1.0,-1.5, -2990, 39.0, 1.0, 0.0, 34.5, 14, 1.0 };
  unsigned array_int2[]={11,1,2,3,4,5,6,7,8,9,1};
  Print(array_double3, 11);
  unsigned index2[]={0,1,2,3,4,5,6,7,8,9,10};
  QuickSort(array_double3, 0,10, index2, array_int2);
  cout<<"***sorted array"<<endl;
  Print(array_double3,11);
  cout<<"***sorted index"<<endl;
  Print(index2,11);
  cout<<"****sorted secondary array"<<endl;
  Print(array_int2, 11);

  cout<<"####now testing reverse array......"<<endl;
  cout<<"\t before.."<<endl;
  Print(array_int,4);
  Reverse(array_int, 4);
  cout<<"\t after...."<<endl;
  Print(array_int, 4);

  //cout<<"#####new test of reverse:"<<endl;
  Print(array_double, 4);
  Reverse(array_double, 4);
  Print(array_double, 4);

  Reverse(vec_int);
  //Print(vec_int);
  Reverse(vec_double);
  //Print(vec_double, 4);


  //testing
  cout<<"%%%%%%testing copy element function"<<endl;
  unsigned* dest=new unsigned[4];
  unsigned index_array[]={2,1,2};
  cout<<"--source:"<<endl;
  Print(array_int,4);
  if(CopyElements(array_int, 4, dest, 4, index_array,3))
    {
      cout<<"--Done"<<endl;
      Print(dest, 4);
    }
  else
    {
      cout<<"ERROR, can not copy"<<endl;
    }

  //testing align with constraints fixed left remove right errors
  cout<<"%%%%%%%%testing alignment function 1 fixed left remove"<<endl;
  double error_cost=5;
  string seq1="CAAGCCCCTTCATCAGGAACTGGTCCGTCGGGTCCCGGCGACACGGGGGTCTCCACGAGAACCTCCTCCCACGATCCCCCTTCTGGCTACCCGGGAACCACCTCCGACTCCTCTGCCACTGGCGCCAGGGAACCGGGGTCTGCAGGTATCGCATCATCGGTTTGGGACCCCATTGGGGAGGTTTTTAGCACTAGGGCCGACTCACCAGAGAGCGTGTCATTATGTGACGGCACAGGAGCCGAGCGGGGTACATGAGACGCAACTATGGTGACTAACGGGATATCAATGAGCAGACAGAGA";
  string seq2=    "GACTTCTCTGCCACTGGTAACAGGGAACCGGGGTCTGTAGTTTTCGTAGT";
  
  /*string seq1="GACTCCTCTGCCACTGGCGCCAGGGAACCGGGGTCTGCAGGTATCGCATCATC";
  string seq2="GACTCCTCTGTCACTGGTCCCACGGTGCCGGGGTCTCTAGCTTCATGGTCATC";
  */
  unsigned maximum_errors=10;
  unsigned n_errors=0;
  unsigned* error_positions=new unsigned [maximum_errors];
  unsigned alength;
  /*  alength=align_with_constraints_fixed_left_remove_right_errors(seq1,seq2, maximum_errors, error_cost, &n_errors, error_positions);

  cout<<"\t aligned length: "<<alength<<";n_error:"<<n_errors<<endl;
  cout<<"\t error positions:";
  for(unsigned i =0;i<maximum_errors;i++)
    {
      cout<<error_positions[i]<<",";
    }
  cout<<endl;
  cout<<"Done"<<endl;*/
  							
  //testing align with constraints fast left
  cout<<"#####################start doing the second testing of alignment functions"<<endl;
  unsigned* align_positions=new unsigned[2];
  alength=align_with_constraints_fast_left(seq1, seq2, maximum_errors, error_cost, align_positions, &n_errors, error_positions);

  cout<<"\t aligned length: "<<alength<<";n_error:"<<n_errors<<endl;
  cout<<"\t error positions:";
  for(unsigned i =0;i<maximum_errors;i++)
    {
      cout<<error_positions[i]<<",";
    }
  
  cout<<endl;

  cout<<"\talign positions:"<<align_positions[0]<<","<<align_positions[1]<<endl;
  cout<<"Done"<<endl;

  //now testing Unique
  cout<<"testing unique function:"<<endl;
  unsigned test_array1[]={1,2,3,3,3,3,3,55,55,7,6,2,1};
  //intialize the output array
  unsigned size_test_array=13;
  unsigned* test_out=new unsigned [size_test_array];
  unsigned* test_out_index=new unsigned[size_test_array];
  for(unsigned i=0;i<size_test_array;i++)
    {
      test_out_index[i]=i;
    }
  unsigned size_unique;
  Unique(test_array1, size_test_array, test_out, test_out_index, size_unique);
  cout<<"\treturn "<<size_unique<<" elements."<<endl;
  
  for(unsigned i=0;i<size_test_array;i++)
    {
      cout<<test_out[i]<<"-"<<test_out_index[i]<<",";
    }
  cout<<endl;
  //clear the memory
  delete[] error_positions;
  									 
  cout<<"!!!!!DONE!!!!!!"<<endl;
  
  return 0;
}

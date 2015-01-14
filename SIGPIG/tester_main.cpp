#include <iostream>
#include <vector>
#include <cstring>

#include "MatrixFunctions.hpp"
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
  cout<<"\tmedian:first:"<<0<<";last:"<<7<<";middle:"<<3<<";median:"<<GetMedianIndex(array_double2,0,7,3)<<endl;
  cout<<"\t\tpivot value:"<<array_double2[GetMedianIndex(array_double2,0,7,3)]<<endl;

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

  
  Reverse(array_int, 4);
  Print(array_int, 4);
  Reverse(array_double, 4);
  Print(array_double, 4);

  Reverse(vec_int);
  //Print(vec_int);
  Reverse(vec_double);
  //Print(vec_double, 4);
  
  return 0;
}

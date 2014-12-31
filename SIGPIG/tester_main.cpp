#include <iostream>
#include <vector>

#include "MatrixFunctions.hpp"
//this is a testing module for testing various functions

using namespace std;
int main()
{
  //first testing matrixFunctions

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
  

  cout<<"return max : "<<max_mf(array_int, 4)<<endl;
  cout<<"return max : "<<max_mf(array_double, 4)<<endl;
  
  cout<<"return max : "<<max_mf(vec_int)<<endl;
  cout<<"return max : "<<max_mf(vec_double)<<endl;

  cout<<"return min : "<<min_mf(array_int, 4)<<endl;
  cout<<"return min : "<<min_mf(array_double, 4)<<endl;
  //const vector<unsigned> constVec=vec_int;
  cout<<"return min : "<<min_mf(vec_int)<<endl;
  cout<<"return min : "<<min_mf(vec_double)<<endl;
  return 0;
}

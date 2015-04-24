#include <iostream>

#include "Poisson.hpp"

using namespace std;
int main()
{

  cout<<"****Testing Poisson distribution******"<<endl;

  Poisson pois;

  //now calling the function
  double p=pois.Density(10, 500);
  cout<<"calling pois.Density(10,11):"<<p<<endl;
  cout<<"***DONE*****"<<endl;
  return 0;
}

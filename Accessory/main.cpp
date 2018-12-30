#include <iostream>
#include <string.h>
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

  char c[1000];
  cout<<"the string is: "<<c<<endl;
  cout<<"the len is :"<<strlen(c)<<endl;
  cout<<"the char is:"<<(int)(c[0])<<endl;
  
  return 0;
}

#include <iostream>
#include <string.h>
#include "Poisson.hpp"
#include <stdio.h>
#include <stdlib.h>

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

  cout<<"Testing memcpy........."<<endl;
  unsigned dst[]={1,2,3,4,5,7,8,9,10};
  cout<<dst[1]<<endl;
  cout<<"Copy/shift 2-3 to 3-4......"<<endl;
  memcpy(&(dst[1]), &(dst[2]), sizeof(unsigned)*2);
	
  cout<<"[0]:"<<dst[0]<<";[1]:"<<dst[1]<<";[2]:"<<dst[2]<<";[3]:"<<dst[3]<<";[4]:"<<dst[4]<<endl;
  
  return 0;
}

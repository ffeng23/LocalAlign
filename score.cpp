#include <stdlib.h>
#include "score.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;
//all matrix data are from ncbi ftp://ftp.ncbi.nih.gov/blast/matrices/
//so far supports only BLOSUM50 and nuc4.4

//********************NUC44***********************
//the actually nuc44 matrix
//the scale information is from matlab. not sure where this is coming from, but
//the score are from ncbi.
//order is A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
static int nuc44Int_1d[] ={ 5,  -4,  -4,  -4,  -4,   1,   1,  -4,  -4,   1,  -4,  -1,  -1,  -1,  -2,
			    -4,  5,  -4,  -4,  -4,   1,  -4,   1,   1,  -4,  -1,  -4,  -1,  -1,  -2,
			    -4, -4,   5,  -4,   1,  -4,   1,  -4,   1,  -4,  -1,  -1,  -4,  -1,  -2,
			    -4, -4,  -4,   5,   1,  -4,  -4,   1,  -4,   1,  -1,  -1,  -1,  -4,  -2,
			    -4, -4,   1,   1,  -1,  -4,  -2,  -2,  -2,  -2,  -1,  -1,  -3,  -3,  -1,
			    1,   1,  -4,  -4,  -4,  -1,  -2,  -2,  -2,  -2,  -3,  -3,  -1,  -1,  -1,
			    1,  -4,   1,  -4,  -2,  -2,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -1,  -1,
			    -4,   1,  -4,   1,  -2,  -2,  -4,  -1,  -2,  -2,  -1,  -3,  -1,  -3,  -1,
			    -4,   1,   1,  -4,  -2,  -2,  -2,  -2,  -1,  -4,  -1,  -3,  -3,  -1,  -1,
			    1,  -4,  -4,   1,  -2,  -2,  -2,  -2,  -4,  -1,  -3,  -1,  -1,  -3,  -1,
			    -4,  -1,  -1,  -1,  -1,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -2,  -2,  -1,
			    -1,  -4,  -1,  -1,  -1,  -3,  -1,  -3,  -3,  -1,  -2,  -1,  -2,  -2,  -1,
			    -1,  -1,  -4,  -1,  -3,  -1,  -3,  -1,  -3,  -1,  -2,  -2,  -1,  -2,  -1 , 
			    -1,  -1,  -1,  -4,  -3,  -1,  -1,  -3,  -1,  -3,  -2,  -2,  -2,  -1,  -1,
			    -2,  -2,  -2,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1
		            };
static const int* nuc44Int[] ={nuc44Int_1d,nuc44Int_1d+15,nuc44Int_1d+30,nuc44Int_1d+45,nuc44Int_1d+60,nuc44Int_1d+75,
			 nuc44Int_1d+90,nuc44Int_1d+105, nuc44Int_1d+120, nuc44Int_1d+135,nuc44Int_1d+150, 
			 nuc44Int_1d+165, nuc44Int_1d+15*12,nuc44Int_1d+15*13, nuc44Int_1d+15*14 
                };
const char ScoreMatrix::NucAlphabet[]={ 'A' ,  'T' ,  'G'  ,'C',   'S',   'W',   'R' ,  'Y',   'K',   'M',   'B' ,  'V' ,  'H' , 'D',   'N'};
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix nuc44(nuc44Int,0.277316, ScoreMatrix::NucAlphabet,15);//again not sure where this scale come from by Matlab.
//****************************************

//*********************************************
//BLOSUM50, aa
//the scale information is from ncbi. not sure where this is coming from, but
//the score are from ncbi.
//order is A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X,*
static int blosum50_1d[] ={  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0,-2,-1,-1,-5,
			    -2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3,-1, 0,-1,-5,
			    -1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3, 4, 0,-1,-5,
			    -2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4, 5, 1,-1,-5,
			    -1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-3,-3,-2,-5,
			     -1,1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3, 0, 4,-1,-5,
			     -1,0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3, 1, 5,-1,-5,
			      0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4,-1,-2,-2,-5,
			     -2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4,0,0,-1,-5,
			     -1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4,-4,-3,-1,-5,
			     -2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1,-4,-3,-1,-5,
			     -1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3,0,1,-1,-5,
			     -1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1,-3,-1,-1,-5,
			     -3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1,-4,-4,-2,-5,
			     -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3,-2,-1,-2,-5,
			     1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2,0,0,-1,-5,
			     0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0,0,-1,0,-5,
			     -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,2,-3,-5,-2,-3,-5,
			     -2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1,-3,-2,-1,-5,
			     0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5,-4,-3,-1,-5,
			     -2,-1,4,5,-3,0,1,-1,0,-4,-4,0,-3,-4,-2,0,0,-5,-3,-4,5,2,-1,-5,
			     -1,0,0,1,-3,4,5,-2,0,-3,-3,1,-1,-4,-1,0,-1,-2,-2,-3,2,5,-1,-5,
			     -1,-1,-1,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1,0,-3,-1,-1,-1,-1,-1,-5,
			     -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,1 };
static const int* blosum50Int[] ={blosum50_1d,blosum50_1d+24,blosum50_1d+24*2,blosum50_1d+24*3, blosum50_1d+24*4,
			    blosum50_1d+24*5,blosum50_1d+24*6,blosum50_1d+24*7,blosum50_1d+24*8,blosum50_1d+24*9,
			    blosum50_1d+24*10,blosum50_1d+24*11,blosum50_1d+24*12,blosum50_1d+24*13,blosum50_1d+24*14,
			    blosum50_1d+24*15,blosum50_1d+24*16,blosum50_1d+24*17,blosum50_1d+24*18,blosum50_1d+24*19,
			    blosum50_1d+24*20,blosum50_1d+24*21,blosum50_1d+24*22,blosum50_1d+24*23};



ScoreMatrix blosum50(blosum50Int,log(2)/3,ScoreMatrix::AaAlphabet,24);//again not sure where this scale come from, but I got it by google.
                                           //it is from NCBI C toolkit cross reference
                                           //http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM50
const char ScoreMatrix::AaAlphabet[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*'};
//****************************************************************



ScoreMatrix::ScoreMatrix(const int*  matrix[ ] , const double& scale, const char* alphabet, const int& alphabet_array_size)
  :
  _matrix(matrix),_scale(scale),
  _alphabet(alphabet),_scaledScoreFlag(false), c_alphabet_array_size(alphabet_array_size)
{
  //do nothing 
}
const int ScoreMatrix::GetAlphabetLength()
{
  return this->c_alphabet_array_size;
}
const int** ScoreMatrix::GetScoreMatrix()
{
  return this->_matrix;
}
  
double ScoreMatrix::GetScale()
{
  return this->_scale;
}

ScoreMatrix::~ScoreMatrix()
{
  this->_matrix=NULL;
}
const char* ScoreMatrix::GetAlphabet()
{
  return this->_alphabet;
}
  
double ScoreMatrix::GetScore(const char& first, const char& second)
{
  //look through
  
  int firstIndex=-1, secondIndex=-1;
  char firstUp=toupper(first);
  char secondUp=toupper(second);
  for(int i=0;i<c_alphabet_array_size;i++)
    {
      
      if(firstUp==_alphabet[i])
	{
	  firstIndex=i;
	}
      if(secondUp==_alphabet[i])
	{
	  secondIndex=i;
	}
      
      if(firstIndex!=-1&&secondIndex!=-1)
	{
	  break;
	}
    }
 
  if(firstIndex==-1)
	{
	  cout<<"********WARNING:can not find the letter \'"<<first<<"\', using \'"<<_alphabet[c_alphabet_array_size-1]<<"\' instead."<<endl; 
	  firstIndex=c_alphabet_array_size-1;
	}
  if(secondIndex==-1)
	{
	  cout<<"********WARNING:can not find the letter \'"<<first<<"\', using \'"<<_alphabet[c_alphabet_array_size-1]<<"\' instead."<<endl; 
	  secondIndex=c_alphabet_array_size-1;
	}

  if(_scaledScoreFlag)
    {
      //cout<<"in no scaled case:**"<<endl;
      return _matrix[firstIndex][secondIndex]/_scale;
    }
  else
    {
      //cout<<"in yes scaled case:"<<endl;
      return _matrix[firstIndex][secondIndex];
    }
}

void ScoreMatrix::SetScaledScoreFlag(const bool& flag)
{
  this->_scaledScoreFlag=flag;
}

//***************************************************************
//the tsm1 matrix
//order is A   T   C H K
static int tsm1Int_1d[] ={  10,  -4,  -4,  -4,  -4,
			     -4,   5,  -4,  -4,   5,
			     -4,  -4,   5,  -4,   6,
			     -4,  -4,  -4,  -4,  -4,
			     -4,   5,   6,  -4,   5,           };
static const int* tsm1Int[] ={tsm1Int_1d,tsm1Int_1d+5*1,tsm1Int_1d+5*2,tsm1Int_1d+5*3,tsm1Int_1d+5*4                };
static const char tsm1Alphabet[]={ 'A' , 'C',   'H',    'K', 'T'   };
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix tsm1(tsm1Int,1, tsm1Alphabet,5);//again not sure where this scale come from by Matlab.
//***************

//***************************************************
//the tsm2 matrix
//order is A   T   C  G
static int tsm2Int_1d[] ={  10,  -9,  -9,  -9,
			    -9,  10,   -9,   -9,
			    -9,  -9,   10,   -9,
			    -9,  -9,   -9,   10
};
static const int* tsm2Int[] ={tsm2Int_1d,tsm2Int_1d+4*1,tsm2Int_1d+4*2,tsm2Int_1d+4*3                };
static const char tsm2Alphabet[]={ 'A' , 'T', 'C',   'G'   };
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix tsm2(tsm2Int,1, tsm2Alphabet,4);
//this matrix is used to test the example in the reference:Barton,G. 1993. CABIO. An efficient Algorithm to locate all locally ptimal alignments between two sequences alowing for gaps. Vol 9. no.6. 1993 P729-34
//******************************************

//*************************NUC44 high penalty of MisMatch
//the actually nuc44 matrix modified with high penalty of mismatch
//the scale information is kept same as the original one from NUC44
//the score are from ncbi.
//order is A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
static int factor=4;
static int nuc44HPInt_1d[] ={ 5,  -4*factor,-4*factor,  -4*factor,  -4*factor, 1,   1,  -4*factor,  -4*factor,   1,    -4*factor,  -1*factor,  -1*factor,  -1*factor,  -2*factor,
			    -4*factor,  5,  -4*factor,  -4*factor,  -4*factor,   1,  -4*factor,   1,   1,  -4*factor,  -1*factor,  -4*factor,  -1*factor,  -1*factor,  -2*factor,
			    -4*factor, -4*factor,   5,  -4*factor,   1,  -4*factor,   1,  -4*factor,   1,  -4*factor,  -1*factor,  -1*factor,  -4*factor,  -1*factor,  -2*factor,
			    -4*factor, -4*factor,  -4*factor,   5,   1,  -4*factor,  -4*factor,   1,  -4*factor,   1,  -1*factor,  -1*factor,  -1*factor,  -4*factor,  -2*factor,
			    -4*factor, -4*factor,   1,   1,  -1*factor,  -4*factor,  -2*factor,  -2*factor, -2*factor,  -2*factor,  -1*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,
			    1,   1,  -4*factor,  -4*factor,  -4*factor,  -1*factor,  -2*factor,  -2*factor,  -2*factor,  -2*factor,  -3*factor,  -3*factor,  -1*factor,  -1*factor,  -1*factor,
			    1,  -4*factor,   1,  -4*factor,  -2*factor,  -2*factor,  -1*factor,  -4*factor,  -2*factor,  -2*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,  -1*factor,
			    -4*factor,   1,  -4*factor,   1,  -2*factor,  -2*factor,  -4*factor,  -1*factor,  -2*factor,  -2*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,
			    -4*factor,   1,   1,  -4*factor,  -2*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,  -4*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,  -1*factor,
			    1,  -4*factor,  -4*factor,   1,  -2*factor,  -2*factor,  -2*factor,  -2*factor,  -4*factor,  -1*factor,  -3*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,
			    -4*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,
			    -1*factor,  -4*factor,  -1*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,  -2*factor,  -1*factor,  -2*factor,  -2*factor,  -1*factor,
			    -1*factor,  -1*factor,  -4*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,  -2*factor,  -2*factor,  -1*factor,  -2*factor,  -1*factor , 
			    -1*factor,  -1*factor,  -1*factor,  -4*factor,  -3*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,  -1*factor,
			    -2*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor
		            };
static const int* nuc44HPInt[] ={nuc44HPInt_1d,nuc44HPInt_1d+15,nuc44HPInt_1d+30,nuc44HPInt_1d+45,nuc44HPInt_1d+60,nuc44HPInt_1d+75,
			 nuc44HPInt_1d+90,nuc44HPInt_1d+105, nuc44HPInt_1d+120, nuc44HPInt_1d+135,nuc44HPInt_1d+150, 
			 nuc44HPInt_1d+165, nuc44HPInt_1d+15*12,nuc44HPInt_1d+15*13, nuc44HPInt_1d+15*14 
                };
const char NucHPAlphabet[]={ 'A' ,  'T' ,  'G'  ,'C',   'S',   'W',   'R' ,  'Y',   'K',   'M',   'B' ,  'V' ,  'H' , 'D',   'N'};
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix nuc44HP(nuc44HPInt,0.277316, NucHPAlphabet,15);//again not sure where this scale come from by Matlab.


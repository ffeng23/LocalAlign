#include <stdlib.h>
#include "score.hpp";
#include <Math.h>;


//the actually 4x4 nuc44 matrix
static int nuc44Int [4][4]={{5,-1,-1,-1},
		 {-1,5,-1,-1},
		 {-1,-1,5,-1},
		 {-1,-1,-1,5}
                };

ScoreMatrix nuc44(nuc44Int,math.Log(2));

ScoreMatrix::ScoreMatrix(int*  _matrix[ ] , const double& _scale)
{
  this->_matrix=_matrix;
  this->_scale=_scale;

}

int** ScoreMatrix::GetScoreMatrix()
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


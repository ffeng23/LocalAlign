#ifndef SCORE_HPP
#define SCORE_HPP

//here in this file we define the score matrix for the most commonly 
//used ones



class ScoreMatrix
{
  //empty constructor
public :
  //ScoreMatrix(); this is not allowed in case someone call it and we end up with a uninitialized score matrix
  ScoreMatrix(int** _matrix , const double& _scale);

  //destructor
  ~ScoreMatrix();

  int** GetScoreMatrix();  
  double GetScale();

  //member
private:		       
  int** _matrix;
  double _scale;

};

extern ScoreMatrix nuc44 ;


#endif

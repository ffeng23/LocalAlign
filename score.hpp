#ifndef SCORE_HPP
#define SCORE_HPP

//here in this file we define the score matrix for the most commonly 
//used ones



class ScoreMatrix
{
  //empty constructor
public :
  //ScoreMatrix(); this is not allowed in case someone call it and we end up with a uninitialized score matrix
  ScoreMatrix(const int* matrix[] , const double& scale, const char* alphabet, const int &alphabet_array_size);
  
  //we did not define the copy constructor, since the default one would be good.

  //destructor
  ~ScoreMatrix();

  //getter
  const int** GetScoreMatrix();  
  double GetScale();
  const char* GetAlphabet();
  void SetScaledScoreFlag(const bool& Flag);
  
  double GetScore(const char& first, const char& second);
  const int GetAlphabetLength();
  //static members
  //this is the defualt static alphabet array for usering to use
  static const char AaAlphabet[]; //default amino acid alphabetic array
  static const char NucAlphabet[];//default Nucletide array.

  //member
private:		       
  const int** _matrix;
  double _scale;
  const char* _alphabet;
  bool _scaledScoreFlag;//this one is used to indicated wether we return
  //scaled score or not. default using unscaled score
  const int c_alphabet_array_size;
};

extern ScoreMatrix nuc44 ;

extern ScoreMatrix blosum50;

extern ScoreMatrix tsm1;
extern ScoreMatrix tsm2;
#endif

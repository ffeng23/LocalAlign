#ifndef BASE_MODEL_HPP
#define BASE_MODEL_HPP

#include "Assigns.hpp"
#include "Counter.hpp"
#include "../SIGPIG/AlignmentSettings.hpp"
//this is the class define the abstract base model for fitting

//abstract class -- VDJ base model, define some generic method
class BaseModel
{
public:
  BaseModel();
  virtual ~BaseModel()=0;
  
  void ReadSettings(const AlignmentSettings& _asettings);

  virtual bool ValidateModel()=0;//done know what to validate.
  //check for some "important" model parameters and make sure they are set!!!
  //who know what I will be check.
  //the inherited class will know what to do.
  
  virtual bool Normalize()=0;//normalize the model
  //the inherited class will decide what to do.

  //initialize the input counter, set up its parameter and make it ready to
  //do the model counting.
  //the outer caller need to make the couter available as well as
  //the counter is base class and will be determine by the 
  //inherited call to decide which one to use
  virtual bool InitializeCounter(Counter& _c)=0 const;
  
  //the user will supply the model to be populated
  //polymorphism here!! 
  virtual void GetModelFromCounter(const Counter& _c)=0; 
  
  //initialize the assign and set up the parameter
  //the caller needs to make the assigns available
  //inherited class will decide which one to work on
  //polymorphism
  virtual void InitializeAssign(Assigns& _a)=0 const;
  
  //_c is output
  virtual void UpdateCounter(const Assigns& _a, Counter& _c)=0 const;

  //sum counter
  virtual Counter SumCounter(const Counter& _c1, const Counter& _c2)=0 const;
  
  
};


#endif

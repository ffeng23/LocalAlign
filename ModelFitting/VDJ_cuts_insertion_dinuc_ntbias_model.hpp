#ifndef VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_HPP
#define VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_HPP

#include "BaseModel.hpp"
//inherited class
class VDJ_cuts_insertion_dinuc_ntbias_model:public BaseModel
{
public:
  VDJ_cuts_insertion_dinuc_ntbias_model();

  virtual ~VDJ_cuts_insertion_dinuc_ntbias_model();

  virtual bool ValidateModel();

  //define the members

};


#endif


#include <cstring>
#include <sstream>
#include "../SequenceString.hpp"
#include "genomicSegments.hpp"
#include "Alignment.hpp"
#include "MatrixFunctions.hpp"
#include "../string_ext.hpp"
#include "AlignmentSettings.hpp"

#include "Alignment_D.hpp"


unsigned Alignment_D::allele_order []={ 0,
					  1, 2, 3, 4, 5, 6, 7, 8, 9,10,
					  11,12,13,14,15,16,17,18,19,20,
					  21,22,23,24,25,26,27,28,29,30,
					  31,32,33
					  };

unsigned Alignment_D:n_D_alleles=AlignmentSettings::n_D_alleles;


Alignment_D::Alignment_D(): numOfAligned(NULL), align_length(NULL), score(NULL),
			    n_errors(NULL), error_positions(NULL),
			    excess_error_positions_left(NULL),
			    excess_error_positions_right(NULL),
			    align_position_left(NULL), align_position_right(NULL),
			    deletions_left(NULL), deletions_right(NULL),
			    p_region_max_length_left(NULL),
			    p_region_max_length_right(NULL)
{
  //empty one
}

//copy constructor, we need to do a deep copy
Alignment_D::Alignment_D(const Alignment_D& _aod):
             numOfAligned(NULL), align_length(NULL), score(NULL),
			    n_errors(NULL), error_positions(NULL),
			    excess_error_positions_left(NULL),
			    excess_error_positions_right(NULL),
			    align_position_left(NULL), align_position_right(NULL),
			    deletions_left(NULL), deletions_right(NULL),
			    p_region_max_length_left(NULL),
			    p_region_max_length_right(NULL)
{
  unsigned sizeOfUnsigned=sizeof(unsigned)/sizeof(char);
  if(_aod.numOfAligned!=NULL)
    {
      numOfAligned=new unsigned[_aod.n_D_alleles];
      memcpy(numOfAligned, _aod.numOfAligned, n_D_alleles*sizeOfUnsigned);
    }

  if(_aod.align_length!=NULL)
    {
      align_length=new unsigned* [n_D_alleles];
      //memcpy(align_length, _aod.align_length, n_D_alleles*sizeOfUnsigned);
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  align_length[i]=NULL;
	  if(_aod.align_length[i]!=NULL)
	    {
	      align_length[i]=new unsigned[numOfAligned[i]];
	      memcpy(align_length[i], _aod.align_length[i], numOfAligned[i]*sizeOfUnsigned);
	    }
	}
    }//end of align_length

  //score **
  if(_aod.score!=NULL)
    {
      score=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  score[i]=NULL;
	  if(_aod.score[i]!=NULL)
	    {
	      score[i]=new unsigned [numOfAligned[i]];
	      memcpy(score[i], _aod.score[i], numOfAligned[i]*sizeOfUnsigned);
	    }
	}
    }//end of socre

  //n_errors **
  if(_aod.n_errors!=NULL)
    {
      n_errors=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  n_errors[i]=NULL;
	  if(_aod.n_errors[i]!=NULL)
	    {
	      n_errors[i]=new unsigned [numOfAligned[i]];
	      memcpy(n_errors[i], _aod.n_errors[i], numOfAligned[i]*sizeOfUnsigned);
	    }
	}
    }//end of n_errors

  //error_positions ***
  if(_aod.error_positions!=NULL)
    {
      error_positions=new unsigned**[n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  error_positions[i]=NULL;
	  if(_aod.error_positions[i]!=NULL)
	    {
	      error_positions[i]=new unsigned*[numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  error_positions[i][j]=NULL;
		  
		  if(_aod.error_positions[i][j]!=NULL)
		    {
		      error_positions[i][j]=new unsigned[n_errors[i][j]];
		      memcpy(error_positions[i][j], _aod.error_positions[i][j],
			     n_errors[i][j]*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of error_positions

  //excess_error_postions_left **
  if(_aod.excess_error_positions_left!=NULL)
    {
      excess_error_positions_left=new unsigned**[n_D_alleles];
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  excess_error_positions_left[i]=NULL;
	  if(_aod.excess_error_positions_left[i]!=NULL)
	    {
	      excess_error_positions_left[i]=new unsigned* [numberOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  excess_error_positions_left[i][j]=NULL;
		  if(_aod.excess_error_positions_left[i][j]!=NULL)
		    {
		      excess_error_positions_left[i][j]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
		      memcpy(excess_error_positions_left[i][j], _aod.excess_error_positions_left[i][j], AlignmentSettings::negative_excess_deletions_max * sizeOfUnsigned);
		    }
		  
		}
	    }
	}
    }//excess_error_positions_left

  //excess_error_postions_right **
  if(_aod.excess_error_positions_right!=NULL)
    {
      excess_error_positions_right=new unsigned**[n_D_alleles];
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  excess_error_positions_right[i]=NULL;
	  if(_aod.excess_error_positions_right[i]!=NULL)
	    {
	      excess_error_positions_right[i]=new unsigned* [numberOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  excess_error_positions_right[i][j]=NULL;
		  if(_aod.excess_error_positions_right[i][j]!=NULL)
		    {
		      excess_error_positions_right[i][j]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
		      memcpy(excess_error_positions_right[i][j], _aod.excess_error_positions_right[i][j], AlignmentSettings::negative_excess_deletions_max * sizeOfUnsigned);
		    }
		  
		}
	    }
	}
    }//excess_error_positions_right

  //align_position_left **
  if(_aod.align_length_left!=NULL)
    {
      align_length_left=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  align_length_left[i]=NULL;
	  if(_aod.align_length_left[i]!=NULL);
	  {
	    align_length_left[i]=nw unsigned [numOfAligned];
	    memcpy(align_length_left[i], _aod.align_length_left[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of align_length_left**

  //align_position_right **
  if(_aod.align_length_right!=NULL)
    {
      align_length_right=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  align_length_right[i]=NULL;
	  if(_aod.align_length_right[i]!=NULL);
	  {
	    align_length_right[i]=nw unsigned [numOfAligned];
	    memcpy(align_length_right[i], _aod.align_length_right[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of align_length_left**

  //deletion_left **
  if(_aod.deletions_left!=NULL)
    {
      deletions_left=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  deletions_left[i]=NULL;
	  if(_aod.deletions_left[i]!=NULL);
	  {
	    deletions_left[i]=nw unsigned [numOfAligned];
	    memcpy(deletions_left[i], _aod.deletions_left[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of deletions_left**

  //deletion_left **
  if(_aod.deletions_right!=NULL)
    {
      deletions_right=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  deletions_right[i]=NULL;
	  if(_aod.deletions_right[i]!=NULL);
	  {
	    deletions_right[i]=nw unsigned [numOfAligned];
	    memcpy(deletions_right[i], _aod.deletions_right[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of deletions_right**

  //p_region_max_length_left
  if(_aod.p_region_max_length_left!=NULL)
    {
      p_region_max_length_left=new unsigned**[n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  p_region_max_length_left[i]=NULL;
	  if(_aod.p_region_max_length_left[i]!=NULL)
	    {
	      p_region_max_length_left[i]=new unsigned* [numOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  p_region_max_length_left[i][j]=NULL;
		  if(_aod.p_region_max_length_left[i][j]!=NULL)
		    {
		      p_region_max_length_left[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletions];
		      memcpy(p_region_max_length_left[i][j], _aod.p_region_max_length_left[i][j], (1+AlignmentSettings::D_maximum_deletions)*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of p_region_max_length_left

  //p_region_max_length_left
  if(_aod.p_region_max_length_right!=NULL)
    {
      p_region_max_length_right=new unsigned**[n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  p_region_max_length_right[i]=NULL;
	  if(_aod.p_region_max_length_right[i]!=NULL)
	    {
	      p_region_max_length_right[i]=new unsigned* [numOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  p_region_max_length_right[i][j]=NULL;
		  if(_aod.p_region_max_length_right[i][j]!=NULL)
		    {
		      p_region_max_length_right[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletions];
		      memcpy(p_region_max_length_right[i][j], _aod.p_region_max_length_right[i][j], (1+AlignmentSettings::D_maximum_deletions)*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of p_region_max_length_right

  //Gosh finally done.
}

Alignment_D& Alignment_D::operator = (const Alignment_D& _aod):
  numOfAligned(NULL), align_length(NULL), score(NULL),
			    n_errors(NULL), error_positions(NULL),
			    excess_error_positions_left(NULL),
			    excess_error_positions_right(NULL),
			    align_position_left(NULL), align_position_right(NULL),
			    deletions_left(NULL), deletions_right(NULL),
			    p_region_max_length_left(NULL),
			    p_region_max_length_right(NULL)
{
  if(this==&_aod)
    {
      return *this;
    }

  //start copying
  unsigned sizeOfUnsigned=sizeof(unsigned)/sizeof(char);
  if(_aod.numOfAligned!=NULL)
    {
      numOfAligned=new unsigned[_aod.n_D_alleles];
      memcpy(numOfAligned, _aod.numOfAligned, n_D_alleles*sizeOfUnsigned);
    }

  if(_aod.align_length!=NULL)
    {
      align_length=new unsigned* [n_D_alleles];
      //memcpy(align_length, _aod.align_length, n_D_alleles*sizeOfUnsigned);
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  align_length[i]=NULL;
	  if(_aod.align_length[i]!=NULL)
	    {
	      align_length[i]=new unsigned[numOfAligned[i]];
	      memcpy(align_length[i], _aod.align_length[i], numOfAligned[i]*sizeOfUnsigned);
	    }
	}
    }//end of align_length

  //score **
  if(_aod.score!=NULL)
    {
      score=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  score[i]=NULL;
	  if(_aod.score[i]!=NULL)
	    {
	      score[i]=new unsigned [numOfAligned[i]];
	      memcpy(score[i], _aod.score[i], numOfAligned[i]*sizeOfUnsigned);
	    }
	}
    }//end of socre

  //n_errors **
  if(_aod.n_errors!=NULL)
    {
      n_errors=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  n_errors[i]=NULL;
	  if(_aod.n_errors[i]!=NULL)
	    {
	      n_errors[i]=new unsigned [numOfAligned[i]];
	      memcpy(n_errors[i], _aod.n_errors[i], numOfAligned[i]*sizeOfUnsigned);
	    }
	}
    }//end of n_errors

  //error_positions ***
  if(_aod.error_positions!=NULL)
    {
      error_positions=new unsigned**[n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  error_positions[i]=NULL;
	  if(_aod.error_positions[i]!=NULL)
	    {
	      error_positions[i]=new unsigned*[numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  error_positions[i][j]=NULL;
		  
		  if(_aod.error_positions[i][j]!=NULL)
		    {
		      error_positions[i][j]=new unsigned[n_errors[i][j]];
		      memcpy(error_positions[i][j], _aod.error_positions[i][j],
			     n_errors[i][j]*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of error_positions

  //excess_error_postions_left **
  if(_aod.excess_error_positions_left!=NULL)
    {
      excess_error_positions_left=new unsigned**[n_D_alleles];
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  excess_error_positions_left[i]=NULL;
	  if(_aod.excess_error_positions_left[i]!=NULL)
	    {
	      excess_error_positions_left[i]=new unsigned* [numberOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  excess_error_positions_left[i][j]=NULL;
		  if(_aod.excess_error_positions_left[i][j]!=NULL)
		    {
		      excess_error_positions_left[i][j]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
		      memcpy(excess_error_positions_left[i][j], _aod.excess_error_positions_left[i][j], AlignmentSettings::negative_excess_deletions_max * sizeOfUnsigned);
		    }
		  
		}
	    }
	}
    }//excess_error_positions_left

  //excess_error_postions_right **
  if(_aod.excess_error_positions_right!=NULL)
    {
      excess_error_positions_right=new unsigned**[n_D_alleles];
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  excess_error_positions_right[i]=NULL;
	  if(_aod.excess_error_positions_right[i]!=NULL)
	    {
	      excess_error_positions_right[i]=new unsigned* [numberOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  excess_error_positions_right[i][j]=NULL;
		  if(_aod.excess_error_positions_right[i][j]!=NULL)
		    {
		      excess_error_positions_right[i][j]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
		      memcpy(excess_error_positions_right[i][j], _aod.excess_error_positions_right[i][j], AlignmentSettings::negative_excess_deletions_max * sizeOfUnsigned);
		    }
		  
		}
	    }
	}
    }//excess_error_positions_right

  //align_position_left **
  if(_aod.align_position_left!=NULL)
    {
      align_position_left=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  align_position_left[i]=NULL;
	  if(_aod.align_position_left[i]!=NULL);
	  {
	    align_position_left[i]=nw unsigned [numOfAligned];
	    memcpy(align_position_left[i], _aod.align_position_left[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of align_length_left**

  //align_position_right **
  if(_aod.align_position_right!=NULL)
    {
      align_position_right=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  align_position_right[i]=NULL;
	  if(_aod.align_position_right[i]!=NULL);
	  {
	    align_position_right[i]=new unsigned [numOfAligned];
	    memcpy(align_position_right[i], _aod.align_position_right[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of align_length_left**

  //deletion_left **
  if(_aod.deletions_left!=NULL)
    {
      deletions_left=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  deletions_left[i]=NULL;
	  if(_aod.deletions_left[i]!=NULL);
	  {
	    deletions_left[i]=new unsigned [numOfAligned];
	    memcpy(deletions_left[i], _aod.deletions_left[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of deletions_left**

  //deletion_left **
  if(_aod.deletions_right!=NULL)
    {
      deletions_right=new unsigned* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  deletions_right[i]=NULL;
	  if(_aod.deletions_right[i]!=NULL);
	  {
	    deletions_right[i]=new unsigned [numOfAligned];
	    memcpy(deletions_right[i], _aod.deletions_right[i], numOfAligned*sizeOfUnsigned);
	  }
	}
    }//end of deletions_right**

  //p_region_max_length_left
  if(_aod.p_region_max_length_left!=NULL)
    {
      p_region_max_length_left=new unsigned**[n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  p_region_max_length_left[i]=NULL;
	  if(_aod.p_region_max_length_left[i]!=NULL)
	    {
	      p_region_max_length_left[i]=new unsigned* [numOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  p_region_max_length_left[i][j]=NULL;
		  if(_aod.p_region_max_length_left[i][j]!=NULL)
		    {
		      p_region_max_length_left[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletions];
		      memcpy(p_region_max_length_left[i][j], _aod.p_region_max_length_left[i][j], (1+AlignmentSettings::D_maximum_deletions)*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of p_region_max_length_left

  //p_region_max_length_left
  if(_aod.p_region_max_length_right!=NULL)
    {
      p_region_max_length_right=new unsigned**[n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  p_region_max_length_right[i]=NULL;
	  if(_aod.p_region_max_length_right[i]!=NULL)
	    {
	      p_region_max_length_right[i]=new unsigned* [numOfAligned];
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  p_region_max_length_right[i][j]=NULL;
		  if(_aod.p_region_max_length_right[i][j]!=NULL)
		    {
		      p_region_max_length_right[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletions];
		      memcpy(p_region_max_length_right[i][j], _aod.p_region_max_length_right[i][j], (1+AlignmentSettings::D_maximum_deletions)*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of p_region_max_length_right

  //Gosh finally done.
  return *this;
}

Alignment_D::~Alignment_D()
{
  //start deleting
  
  if(numOfAligned!=NULL)
    {
      delete[] numOfAligned;
    }

  if(align_length!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{	  
	  if(align_length[i]!=NULL)
	    {
	      delete [] align_length[i];
	    }
	}
      delete [] align_length;
    }//end of align_length

  //score **
  if(score!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(score[i]!=NULL)
	    {
	      delete [] score[i];	      
	    }
	}
      delete [] score;
    }//end of socre

  //n_errors **
  if(n_errors!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(n_errors[i]!=NULL)
	    {
	      delete [] n_errors[i];
	    }
	}
      delete [] n_errors;
    }//end of n_errors

  //error_positions ***
  if(error_positions!=NULL)
    {
     
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  
	  if(error_positions[i]!=NULL)
	    {
	      
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  		  
		  if(error_positions[i][j]!=NULL)
		    {
		      delete [] error_positions[i][j];
		    }
		}
	      delete [] error_positions[i]
	    }
	}
      delete[] error_positions;
    }//end of error_positions

  //excess_error_postions_left **
  if(excess_error_positions_left!=NULL)
    {      
      for(unsigned int i=0;i<n_D_alleles;i++)
	{	  
	  if(excess_error_positions_left[i]!=NULL)
	    {	      
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  if(excess_error_positions_left[i][j]!=NULL)
		    {
		      delete[] excess_error_positions_left[i][j];
		    }
		  
		}
	      delete [] excess_error_positions_left[i];
	    }
	}
      delete[] excess_error_positions_left;
    }//excess_error_positions_left

  //excess_error_postions_right **
  if(excess_error_positions_right!=NULL)
    {
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  if(excess_error_positions_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  if(excess_error_positions_right[i][j]!=NULL)
		    {
		      delete[] excess_error_positions_right[i][j];
		    }		  
		}
	      delete [] excess_error_positions_right[i];
	    }
	}
      delete [] excess_error_positions_right;
    }//excess_error_positions_right

  //align_position_left **
  if(align_position_left!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(align_position_left[i]!=NULL);
	  {
	    delete[] align_position_left[i];
	  }
	}
      delete[] align_position_left;
    }//end of align_length_left**

  //align_position_right **
  if(align_position_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(_aod.align_position_right[i]!=NULL);
	  {
	    delete [] align_position_right[i];
	  }
	}
      delete [] align_position_right;
    }//end of align_length_left**

  //deletion_left **
  if(deletions_left!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(deletions_left[i]!=NULL);
	  {
	    delete [] deletions_left[i];
	  }
	}
      delete [] deletions_left;
    }//end of deletions_left**

  //deletion_left **
  if(deletions_right!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(deletions_right[i]!=NULL);
	  {
	    delete [] deletions_right[i];
	  }
	}
      delete [] deletions_right;
    }//end of deletions_right**

  //p_region_max_length_left
  if(p_region_max_length_left!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(p_region_max_length_left[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  if(p_region_max_length_left[i][j]!=NULL)
		    {
		      delete [] p_region_max_length_left[i][j];
		    }
		}
	      delete [] p_region_max_length_left[i];
	    }
	}
      delete [] p_region_max_length_left;
    }//end of p_region_max_length_left

  //p_region_max_length_left
  if(p_region_max_length_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(p_region_max_length_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned;j++)
		{
		  if(p_region_max_length_right[i][j]!=NULL)
		    {
		      delete[] p_region_max_length_right[i][j];
		    }
		}
	      delete [] p_region_max_length_right[i];
	    }
	}
      delete [] p_region_max_length_right;
    }//end of p_region_max_length_right

  //Gosh finally done.
}

string Alignment_D::toString()
{
  //to be implemented 
}

bool match_D(const SequenceString& _seq,
	     const GenomicD* _genDs, const unsigned& _numOfDSegs,
	     const unsigned& _V_end, const unsigned& _J_start,
	     const unsigned& _flank_length,
	     const ScoreMatrix* _nuc44,
	     const unsigned& _D_minimum_alignment_length,
	     const unsigned& _D_maximum_deletion, 
	     const unsigned& _negative_excess_deletion_max,
	     const unsigned& _D_allowed_errors, 
	     const unsigned& _error_cost,
	     /*output*/ Alignment_D& _D)
{
  unsigned l_seq=_seq.GetLength();
  unsigned start_index=_V_end+1-_flank_length; //plus one because, we started from one past where the V aligned ends;
  if(start_index<0)
    {
      start_index=0;
    }
  unsigned end_index=J_start-1+_flank_length;
  if(end_index>=l_seq)
    {
      end_index=l_seq-1;
    }
  string seqNDN=_seq.GetSequence().substr(start_index, end_index-start_index+1);
  unsigned flank_offset=start_index-1;
}

void DeterminePalindromAndExcessError_D(
					)
{
}

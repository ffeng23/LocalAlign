#include <cstring>
#include <sstream>
#include "../SequenceString.hpp"
#include "genomicSegments.hpp"
#include "Alignment.hpp"
#include "MatrixFunctions.hpp"
#include "../string_ext.hpp"
#include "AlignmentSettings.hpp"
#include "../score.hpp"
#include "../LocalAlignment.hpp"
#include "Alignment_D.hpp"


/*unsigned Alignment_D::allele_order []={ 0,
					  1, 2, 3, 4, 5, 6, 7, 8, 9,10,
					  11,12,13,14,15,16,17,18,19,20,
					  21,22,23,24,25,26,27,28,29,30,
					  31,32,33
					  };
*/
unsigned Alignment_D:n_D_alleles=AlignmentSettings::n_D_alleles;


Alignment_D::Alignment_D(): numOfAligned(NULL), align_length(NULL), score(NULL),
			    n_errors(NULL), error_positions(NULL),
			    excess_error_positions_left(NULL),
			    excess_error_positions_right(NULL),
			    align_position_left(NULL), align_position_right(NULL),
			    deletions_left(NULL), deletions_right(NULL),
			    p_region_max_length_left(NULL),
			    p_region_max_length_right(NULL),
			    allele_order(NULL)
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
	     p_region_max_length_right(NULL),
	     allele_order(NULL)
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
      score=new double* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  score[i]=NULL;
	  if(_aod.score[i]!=NULL)
	    {
	      score[i]=new double [numOfAligned[i]];
	      memcpy(score[i], _aod.score[i], numOfAligned[i]*sizeof(double)/sizeof(char));
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
	    align_length_right[i]=nw unsigned [numOfAligned];
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
  if(_aod.allele_order!=NULL)
    {
      allele.order=new unsigned[n_D_alleles];
      memcpy(allele_order, _aod.allele_order, n_D_alleles*sizeOfUnsigned);
    }
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
  p_region_max_length_right(NULL),allele_order(NUL)
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
      score=new double* [n_D_alleles];
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  score[i]=NULL;
	  if(_aod.score[i]!=NULL)
	    {
	      score[i]=new double [numOfAligned[i]];
	      memcpy(score[i], _aod.score[i], numOfAligned[i]*sizeof(double)/sizeof(char));
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
  //allele_order
  if(_aod.allele_order!=NULL)
    {
      allele.order=new unsigned[n_D_alleles];
      memcpy(allele_order, _aod.allele_order, n_D_alleles*sizeOfUnsigned);
    }
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
  
  //allele_order
  if(allele_order!=NULL)
    delete[] allele_order;
  //Gosh finally done.
}

string Alignment_D::toString()
{
  //to be implemented 
}

//need to implement the error catch system
bool initialize(const unsigned& _n_D_alleles)
{
  n_D_alleles=_n_D_alleles;
  
  //start initialize it
  unsigned sizeOfUnsigned=sizeof(unsigned)/sizeof(char);
  
  numOfAligned=new unsigned[n_D_alleles];
  std::memset(numOfAligned, 0, n_D_alleles*sizeOfUnsigned);
  //set a default value of zero,in case of error. 
  
  align_length=new unsigned* [n_D_alleles];
  
  //score **
  score=new double* [n_D_alleles];
  
  //n_errors **
  n_errors=new unsigned* [n_D_alleles];
  
  //error_positions ***
  error_positions=new unsigned**[n_D_alleles];
  
  //excess_error_postions_left **
  excess_error_positions_left=new unsigned**[n_D_alleles];
  
  //excess_error_postions_right **
  excess_error_positions_right=new unsigned**[n_D_alleles];
  
  //align_position_left **
  align_position_left=new unsigned* [n_D_alleles];
  
  //align_position_right **
  align_position_right=new unsigned* [n_D_alleles];
  
  //deletion_left **
  deletions_left=new unsigned* [n_D_alleles];
  
  //deletion_left **
  deletions_right=new unsigned* [n_D_alleles];
  
  //p_region_max_length_left
  p_region_max_length_left=new unsigned**[n_D_alleles];
  
  //p_region_max_length_left
  p_region_max_length_right=new unsigned**[n_D_alleles];
  
  //allele_order
  allele_order=new unsigned [n_D_alleles];
  return true;
}

bool match_D(const SequenceString& _seq,
	     const GenomicD* _genDs, const unsigned& _numOfDSegs,
	     const unsigned& _V_end, const unsigned& _J_start,
	     const unsigned& _flank_length,
	     const ScoreMatrix* _sm,
	     const unsigned& _D_maximum_deletion, 
	     const unsigned& _negative_excess_deletion_max,
	     const unsigned& _max_aligns,
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
  unsigned D_max_errors=(int)(AlignmentSettings::max_length_D_genes*_sm(1,1)/(_sm(1,1)-_sm(1,2)));  
  //start doing the alignment by calling the localAlign function
  if(seqNDN.size()==0)
    {
      std::memset(_D.numOfAligned, 0, _numOfDSegs,sizeof(unsigned)/sizeof(char));
      return false;
    }
  
  //good, let's do the alignment
  SequenceString ss_seqNDN("temp", seqNDN);
  unsigned curr_numOfAligned;
  _D.initialize(_numOfDSegs);
  AlignmentString* las;
  for(unsigned d=0;d<n_D_alleles;d++)
    {
      LocalAlignment la(&ss_seqNDN, _genDs+d, _nuc44,
			100000/*gapopen*/ ,100000 /*gapextension*/,
			1/*scale*/, _max_aligns);
      //first determine the number of aligned for this current one
      curr_numOfAligned=la.GetNumberOfAlignments();
      las=la.GetAlignmentArr();

      _D.numOfAligned[d]=curr_numOfAligned;
      _D.align_length[d]=new unsigned[curr_numOfAligned];
      _D.score[d]=new double[curr_numOfAligned];
      _D.n_errors[d]=new unsigned[curr_numOfaligned];
      _D.error_positions[d]=new unsigned*[curr_numOfAligned];
      unsigned* temp_error_positions;
      
      _D.aligned_positions_left[d]=new unsigned*[curr_numOfAligned];
      _D.aligned_positions_right[d]=new unsigned*[curr_numOfAligned];
      
      _D.deletions_left[d]=new unsigned* [curr_numOfAligned];
      _D.deletions_right[d]=new unsigned* [curr_numOfAligned];

      for(unsigned i=0;i<curr_numOfAligned;i++)
	{
	  _D.align_length[d][i]=las[i].GetPatternIndexEnd()-las[i].GetPatternIndexStart()+1;
	  _D.score[d][i]=las[i].GetScore();
	  _D.aligned_positions_left[d][i]=las[i].GetPatternIndexStart()+flank_offset+1;
	  _D.aligned_positions_right[d][i]=las[i].GetPatternIndexEnd()+flank_offset+1;
	  _D.deletions_left[d][i]=las[i].GetSubjectIndexStart();
	  _D.deletion_right[d][i]=genDs[d].Get_Seq().GetLength()-las[i].GetSubjectIndexEnd()-1;

	  temp_error_positions=new unsigned[max_n_errors];
	  _D.n_errors[d][i]=findErrors
	    (seqNDN, _genDs+d, _D.aligned_positions_left[d][i],
	     _D.aligned_position_right[d][i],
	     _D.deletions_left[d][i], _D.deletions_right[d][i],
	     max_n_errors, temp_error_positions);
	  //now copy over the elements
	  _D.error_positions[d][i]=new unsigned[_D.n_errors[d][i]];
	  memcpy(_D.error_positions[d][i], temp_error_positions,
		 _D.n_errors[d][i]*sizeof(unsigned)/sizeof(char));
	  
	  //clean up the memory
	  delete [] temp_error_positions;
	}
    }//end of num of D gene segs for loop
  
  //now we are ready with alignment, first sort in order to figure
  //out the D.allele_order, prepare the index array first
  unsigned highestScore=new unsigned[n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      _D.allele_order[i]=i;
      highestScore[i]=_D.score[i][0];
    }
  
  QuickSort(highestScore, 0, _n_D_alleles-1, _D.allele_order,NULL);
  
  //now we are ready to take care of p_nucleotides and negative excess error
  //first need to initialize the p_region array to all zeros
    
}
/*find the error positions for two sequence aligned.
 * the output array has to be allocated by the caller outside 
 */
unsigned findErrors(const string& _seq1, const string& _seq2,
	   const unsigned& pos_start1, const unsigned& pos_end1,
		    const unsigned& post_start2, const unsigned& pos_end2,
		    const unsigned& max_n_errors,
		    /*output*/unsigned* error_position
		    )
{
  //first, do some checking
  if(pos_end1-pos_start1!=pos_end2-pos_start2)
    {
      throw "bad array index, the seq to be check for match/mismatch are not of same length!";
    }
  
  unsigned n_error=0;
  for(unsigned i=0;i<pos_end2-pos_start1+1;i++)
    {
      if(_seq1.at(pos_start1+i)!=_seq2.at(pos_start2+i))
	{	  
	  if(n_error<max_n_errors)
	    error_position[n_error]=i;
	  n_error++;
	}
    }
  
  return n_error;
}

/*Again in this function, the output has been correctly initialized
 *correctly in the caller outside.
 */
void DeterminePalindromAndExcessError_D
( const SequenceString& _seq, const GenomicJ* _genJs,
  /*const unsigned* _ok_order,*/ unsigned** _deletions_left,
  unsigned** _deletions_right,
  const unsigned& _negative_excess_deletions_max, 
  const unsigned& _D_maximum_deletion,
  const unsigned** _align_length, const unsigned& _numOfDSegs,
  const unsigned* _numOfAligned, 
  unsigned** _align_positions_left, unsigned** _align_position_right,
  /*output*/ unsigned*** _p_region_max_length_left, 
  unsigned*** _p_region_max_length_right,
  unsigned*** _excess_error_position_left,
  unsigned*** _excess_error_position_right
  )
{
  bool still_palindrome=true;
  unsigned p=0;
  //now go through to find palindromic nucleotides for each alignment
  for(unsigned d=0;d<_numOfDSegs;d++)
    {
      string target=_genD[d].Get_Sequence();
      //for each aligned in this D seg
      for(unsigned na=0;na<_numOfAligned[d];na++)
	{
	  //getting the left p-nucleotides first
	  //for each possible deletions
	  int nd=_deletions_left[d][na]-_negative_excess_deletions_max;
	  if(nd<0)
	    nd=0;
	  int max_nd=_D_maximum_deletion;
	  if(max_nd>_align_lenth[d][na]+_deletions_left[d][na])
	    max_nd=_align_length[d][na]+_deletions_left[d][na];
	  un
	  for(;nd<=max_nd;nd++)
	    {
	      //for each value of deletions, find longest half-p from the implied end of the gene sequence
	      p=0;
	      still_palindrome=true;
	      //upper bound on p is the length of the match seq or seq to the left of the aligned ones
	      int max_p_length=_align_length[d][na]-(nd-_deletions_left[d][na]);
	      if(max_p_length<0)
		max_p_length=0;
	      if(max_p_length<_align_positions_left[d][na]+(nd-_deletions_left[d][na]))
		{
		  max_p_length=_align_positions_left[d][na]+(nd-_deletions_left[d][na]);
		}
	      while(still_pa&&((signed)p<(signed)max_p_length))
		{
		  still_palindrome=target.at(nd+p)==DnaComplement(_seq.GetSequence().at(_align_positions_left[d][na]-p+nd-_deletions_left[d][na]-1));
		  if(still_palindrome)
		    {
		      p++;
		    }
		}
		
	      _p_region_max_length_left[d][na][1+nd]=p;
	    }//end of for nd<max_nd
	  //start doing excess error 
	  unsigned tempArray[]={_negative_excess_deletions_max, _deletions_left[d][na], align_positions_left[d][na]};
	  unsigned n_excess=min_mf(tempArray,3);
	  unsigned k;
	  unsigned runningIndex_excessError=0;
	  //find the mismatchs
	  for(k=0;k<n_excess;k++)
	    {
	      if(target.at(_deletions_left[d][na]-n_excess+k)!=
		 _seq.GetSequence(_align_position_left[d][na]-n_excess+k))
		{
		  _excess_error_positions_left[d][na][runningIndex_excessError]=
		    _align_positions_left[d][na]-n_excess+k;
		  runningIndex_excessError++;
		}
	    }//end n_excess
	  for(;runningIndex_excessError<_negative_excess_deletions_max;runningIndex_excessError++)
	    {
	      _excess_error_positions_left[d][na][runningIndex_excessError]=0;
	    }

	  //===>>>>doing right side things

	}//end of _numOfAligned
    }//end of outer for loop d : numOfDSegs
}//end of functio find p-nucleotides

#include <cstring>
#include <sstream>
#include <fstream>
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
//unsigned Alignment_D::n_D_alleles=AlignmentSettings::n_D_alleles;
//unsigned Alignnment_D::D_max_errors=

Alignment_D::Alignment_D(): n_D_alleles(0), D_max_errors(0),
			    numOfAligned(NULL), align_length(NULL), score(NULL),
			    n_errors(NULL), error_positions(NULL),			    
			    align_position_left(NULL), align_position_right(NULL),
			    deletions_left(NULL), deletions_right(NULL),
			    p_region_max_length_left(NULL),
			    p_region_max_length_right(NULL),
			    excess_error_positions_left(NULL),
			    excess_error_positions_right(NULL),
			    allele_order(NULL)
{
  //empty one
}

//copy constructor, we need to do a deep copy
Alignment_D::Alignment_D(const Alignment_D& _aod):
             numOfAligned(NULL), align_length(NULL), score(NULL),
	     n_errors(NULL), error_positions(NULL),	     
	     align_position_left(NULL), align_position_right(NULL),
	     deletions_left(NULL), deletions_right(NULL),
	     p_region_max_length_left(NULL),
	     p_region_max_length_right(NULL),
	     excess_error_positions_left(NULL),
	     excess_error_positions_right(NULL),
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
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
	      excess_error_positions_left[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
	      excess_error_positions_right[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
	    align_position_left[i]=new unsigned [numOfAligned[i]];
	    memcpy(align_position_left[i], _aod.align_position_left[i], numOfAligned[i]*sizeOfUnsigned);
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
	    align_position_right[i]=new unsigned [numOfAligned[i]];
	    memcpy(align_position_right[i], _aod.align_position_right[i], numOfAligned[i]*sizeOfUnsigned);
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
	    deletions_left[i]=new unsigned [numOfAligned[i]];
	    memcpy(deletions_left[i], _aod.deletions_left[i], numOfAligned[i]*sizeOfUnsigned);
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
	    deletions_right[i]=new unsigned [numOfAligned[i]];
	    memcpy(deletions_right[i], _aod.deletions_right[i], numOfAligned[i]*sizeOfUnsigned);
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
	      p_region_max_length_left[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  p_region_max_length_left[i][j]=NULL;
		  if(_aod.p_region_max_length_left[i][j]!=NULL)
		    {
		      p_region_max_length_left[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletion];
		      memcpy(p_region_max_length_left[i][j], _aod.p_region_max_length_left[i][j], (1+AlignmentSettings::D_maximum_deletion)*sizeOfUnsigned);
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
	      p_region_max_length_right[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  p_region_max_length_right[i][j]=NULL;
		  if(_aod.p_region_max_length_right[i][j]!=NULL)
		    {
		      p_region_max_length_right[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletion];
		      memcpy(p_region_max_length_right[i][j], _aod.p_region_max_length_right[i][j], (1+AlignmentSettings::D_maximum_deletion)*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of p_region_max_length_right
  if(_aod.allele_order!=NULL)
    {
      allele_order=new unsigned[n_D_alleles];
      memcpy(allele_order, _aod.allele_order, n_D_alleles*sizeOfUnsigned);
    }
  //Gosh finally done.
}

Alignment_D& Alignment_D::operator = (const Alignment_D& _aod)/*:
  numOfAligned(NULL), align_length(NULL), score(NULL),
			    n_errors(NULL), error_positions(NULL),
			    excess_error_positions_left(NULL),
			    excess_error_positions_right(NULL),
			    align_position_left(NULL), align_position_right(NULL),
			    deletions_left(NULL), deletions_right(NULL),
			    p_region_max_length_left(NULL),
			    p_region_max_length_right(NULL),allele_order(NULL)*/
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
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
	      excess_error_positions_left[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
	      excess_error_positions_right[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
	    align_position_left[i]=new unsigned [numOfAligned[i]];
	    memcpy(align_position_left[i], _aod.align_position_left[i], numOfAligned[i]*sizeOfUnsigned);
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
	    align_position_right[i]=new unsigned [numOfAligned[i]];
	    memcpy(align_position_right[i], _aod.align_position_right[i], numOfAligned[i]*sizeOfUnsigned);
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
	    deletions_left[i]=new unsigned [numOfAligned[i]];
	    memcpy(deletions_left[i], _aod.deletions_left[i], numOfAligned[i]*sizeOfUnsigned);
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
	    deletions_right[i]=new unsigned [numOfAligned[i]];
	    memcpy(deletions_right[i], _aod.deletions_right[i], numOfAligned[i]*sizeOfUnsigned);
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
	      p_region_max_length_left[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  p_region_max_length_left[i][j]=NULL;
		  if(_aod.p_region_max_length_left[i][j]!=NULL)
		    {
		      p_region_max_length_left[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletion];
		      memcpy(p_region_max_length_left[i][j], _aod.p_region_max_length_left[i][j], (1+AlignmentSettings::D_maximum_deletion)*sizeOfUnsigned);
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
	      p_region_max_length_right[i]=new unsigned* [numOfAligned[i]];
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  p_region_max_length_right[i][j]=NULL;
		  if(_aod.p_region_max_length_right[i][j]!=NULL)
		    {
		      p_region_max_length_right[i][j]=new unsigned[1+AlignmentSettings::D_maximum_deletion];
		      memcpy(p_region_max_length_right[i][j], _aod.p_region_max_length_right[i][j], (1+AlignmentSettings::D_maximum_deletion)*sizeOfUnsigned);
		    }
		}
	    }
	}
    }//end of p_region_max_length_right
  //allele_order
  if(_aod.allele_order!=NULL)
    {
      allele_order=new unsigned[n_D_alleles];
      memcpy(allele_order, _aod.allele_order, n_D_alleles*sizeOfUnsigned);
    }
  //Gosh finally done.
  return *this;
}

void Alignment_D::ResetData()
{  
  //cout<<"\t****align_length"<<endl;
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
      align_length=NULL;
    }//end of align_length

  //score **
  //cout<<"\t***score"<<endl;
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
      score=NULL;
    }//end of socre
  
  //error_positions ***
  //cout<<"\t****error_positions"<<endl;
  if(error_positions!=NULL)
    {
      //cout<<"\t\t===>not null"<<endl;
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  
	  if(error_positions[i]!=NULL)
	    {
	      //cout<<"\t\t\t==>not null"<<endl;
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  //cout<<"\t\t\t\tn_error["<<i<<"]["<<j<<"]"<<n_errors[i][j]<<endl;	    
		  if(error_positions[i][j]!=NULL)
		    {
		      //		      cout<<"\t\t\t\t%%%%%%%delete it"<<endl;
		      delete [] error_positions[i][j];
		    }
		}
	      delete [] error_positions[i];
	    }
	}
      delete[] error_positions;
      error_positions=NULL;
    }//end of error_positions

  //n_errors **
  //cout<<"\t****n_errors"<<endl;
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
      n_errors=NULL;
    }//end of n_errors

  //excess_error_postions_left **
  //cout<<"\t****excess_error_positions_left"<<endl;
  if(excess_error_positions_left!=NULL)
    {      
      for(unsigned int i=0;i<n_D_alleles;i++)
	{	  
	  if(excess_error_positions_left[i]!=NULL)
	    {	      
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
      excess_error_positions_left=NULL;
    }//excess_error_positions_left

  //excess_error_postions_right **
  //cout<<"\t****excess_error_positions_right"<<endl;
  if(excess_error_positions_right!=NULL)
    {
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  if(excess_error_positions_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
      excess_error_positions_right=NULL;
    }//excess_error_positions_right

  //align_position_left **
  //cout<<"\t****align_position_left"<<endl;
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
      align_position_left=NULL;
    }//end of align_length_left**

  //align_position_right **
  //cout<<"\t****align_position_right"<<endl;
  if(align_position_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(align_position_right[i]!=NULL);
	  {
	    delete [] align_position_right[i];
	  }
	}
      delete [] align_position_right;
      align_position_right=NULL;
    }//end of align_length_left**

  //deletion_left **
  //cout<<"\t****deletions_left"<<endl;
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
      deletions_left=NULL;
    }//end of deletions_left**

  //deletion_left **
  //cout<<"\t****deletion_right"<<endl;
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
      deletions_right=NULL;
    }//end of deletions_right**

  //p_region_max_length_left
  //cout<<"\t****p_left"<<endl;
  if(p_region_max_length_left!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(p_region_max_length_left[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
      p_region_max_length_left=NULL;
    }//end of p_region_max_length_left

  //p_region_max_length_left
  //cout<<"\t****p_right"<<endl;
  if(p_region_max_length_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  if(p_region_max_length_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
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
      p_region_max_length_right=NULL;
    }//end of p_region_max_length_right
  
  //allele_order
  //cout<<"\t****allele_order"<<endl;
  if(allele_order!=NULL)
    {
      delete[] allele_order;
      allele_order=NULL;
    }
  //start deleting
  //cout<<"\t***numOfAligned"<<endl;
  if(numOfAligned!=NULL)
    {
      delete[] numOfAligned;
      numOfAligned=NULL;
    }
  n_D_alleles=0;
  D_max_errors=0;
  //Gosh finally done.
}

Alignment_D::~Alignment_D()
{
  //cout<<"destructing the d alignment object"<<endl;
  ResetData();
  //cout<<"done!!!"<<endl;
}
//for printing/debugging purpose
string Alignment_D::toString() const
{
  stringstream ss;
  //unsigned sizeOfUnsigned=sizeof(unsigned)/sizeof(char);
  ss<<"Alignment D objects:\n";
  ss<<"\tNumber of Aligned:";
  if(numOfAligned!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<numOfAligned[i]<<",";
	}
    }
  ss<<"\n";
  ss<<"\tAligned Length by allele:\n";
  if(align_length!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tallele["<<i<<"]:";
	  if(align_length[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<align_length[i][j]<<",";
		}
	      ss<<"\n";
	    }
	}
    }//end of align_length
  ss<<"\n";

  //score **
  ss<<"\tScore by Allele:\n";
  if(score!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tallele ["<<i<<"]:";
	  if(score[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<score[i][j]<<",";
		}
	      ss<<"\n";
	    }
	}
    }//end of socre
  ss<<"\n";

  //n_errors **
  ss<<"\tnumber of error by allele:\n";
  if(n_errors!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tallele["<<i<<"]:";
	  if(n_errors[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<n_errors[i][j]<<",";
		}
	      ss<<"\n";
	    }
	}
    }//end of n_errors
  ss<<"\n";

  //error_positions ***
  ss<<"\tError Position(# of aligned followed by {error position}):\n";
  if(error_positions!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tAllele["<<i<<"]:";
	  if(error_positions[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  if(n_errors[i][j]==0)
		    continue;
		  
		  ss<<j;
		  if(error_positions[i][j]!=NULL)
		    {
		      ss<<"{";
		      for(unsigned k=0;k<n_errors[i][j];k++)
			{
			  ss<<error_positions[i][j][k]<<",";
			}
		      ss<<"};";
		    }
		  //ss<<"\n";
		}
	    }
	  ss<<"\n";
	}
    }//end of error_positions
  ss<<"\n";

  //excess_error_postions_left **
  ss<<"\tExcess Error Position Left(# of aligned followed by {excess error position}):\n";
  if(excess_error_positions_left!=NULL)
    {
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tAllele["<<i<<"]:";
	  if(excess_error_positions_left[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<j;
		  if(excess_error_positions_left[i][j]!=NULL)
		    {
		      ss<<"{";
		      for(unsigned k=0;k<AlignmentSettings::negative_excess_deletions_max;k++)
			{
			  ss<<excess_error_positions_left[i][j][k]<<",";
			}
		      ss<<"},";
		    }
		  //ss<<"\n";
		}
	    }
	  ss<<"\n";
	}
    }//excess_error_positions_left
  ss<<"\n";
  
  //excess_error_postions_right **
  ss<<"\tExcess Error Position Right (# of aligned followed by {excess error position}):\n";
  if(excess_error_positions_right!=NULL)
    {      
      for(unsigned int i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD Allele["<<i<<"]:";
	  if(excess_error_positions_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<j;
		  if(excess_error_positions_right[i][j]!=NULL)
		    {
		      ss<<"{";
		      for(unsigned k=0;k<AlignmentSettings::negative_excess_deletions_max;k++)
			{
			  ss<<excess_error_positions_right[i][j][k]<<",";
			}
		      ss<<"},";
		    }
		  //ss<<"\n";
		}
	    }
	  ss<<"\n";
	}
    }//excess_error_positions_right
  ss<<"\n";

  //align_position_left **
  ss<<"\tAlign Position left:\n";
  if(align_position_left!=NULL)
    {     
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD allele["<<i<<"]:";
	  if(align_position_left[i]!=NULL);
	  {
	    for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<align_position_left[i][j]<<",";		  
		}
	  }
	  ss<<"\n";	  
	}
    }//end of align_length_left**
  ss<<"\n";
  
  //align_position_right **
  ss<<"\tAlign Position right:\n";
  if(align_position_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD allele["<<i<<"]:";
	  if(align_position_right[i]!=NULL);
	  {
	    for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<align_position_right[i][j]<<",";		  
		}
	  }
	  ss<<"\n";
	}
    }//end of align_length_right**
  ss<<"\n";

  //deletion_left **
  ss<<"\tAlign Deletion Left:\n";
  if(deletions_left!=NULL)
    {      
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD allele["<<i<<"]:";
	  if(deletions_left[i]!=NULL);
	  {
	    for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<deletions_left[i][j]<<",";		  
		}
	  }
	  ss<<"\n";
	}
    }//end of deletions_left**
  ss<<"\n";

  //deletion_right **
  ss<<"\tAlign Deletion Left:\n";
  if(deletions_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD allele["<<i<<"]:";
	  if(deletions_right[i]!=NULL);
	  {
	    for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<deletions_right[i][j]<<",";		  
		}
	  }
	  ss<<"\n";
	}
    }//end of deletions_right**
  ss<<"\n";

  //p_region_max_length_left
  ss<<"\tp_region_max_length_left(# of aligned followed by {p_reg max length}):\n";
  if(p_region_max_length_left!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD allele["<<i<<"]:\n";
	  if(p_region_max_length_left[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<j;
		  if(p_region_max_length_left[i][j]!=NULL)
		    {
		      ss<<"{";
		      for(unsigned k=0;k<1+AlignmentSettings::D_maximum_deletion;k++)
			{
			  ss<<p_region_max_length_left[i][j][k]<<",";
			}
		      ss<<"};";
		    }
		  //ss<<"\n";
		}
	    }
	  ss<<"\n";
	}
    }//end of p_region_max_length_left
  ss<<"\n";

  //p_region_max_length_right
  ss<<"\tp_region_max_length_right(# of aligned followed by {p_reg max length}):\n";
  if(p_region_max_length_right!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<"\t\tD allele["<<i<<"]:\n";
	  if(p_region_max_length_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<numOfAligned[i];j++)
		{
		  ss<<j;
		  if(p_region_max_length_right[i][j]!=NULL)
		    {
		      ss<<"{";
		      for(unsigned k=0;k<1+AlignmentSettings::D_maximum_deletion;k++)
			{
			  ss<<p_region_max_length_right[i][j][k]<<",";
			}
		      ss<<"};";
		    }
		  //ss<<"\n";
		}
	    }
	  ss<<"\n";
	}
    }//end of p_region_max_length_right
  ss<<"\n";
  
  ss<<"\tAllele order:";
  if(allele_order!=NULL)
    {
      for(unsigned i=0;i<n_D_alleles;i++)
	{
	  ss<<allele_order[i]<<",";
	}
      ss<<"\n";
    }

  //Gosh finally done.
  return ss.str();
}

//need to implement the error catch system
bool Alignment_D::initialize(const unsigned& _n_D_alleles)
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
  //cout<<"*&***&&&&&&&:inside matchD"<<endl;
  unsigned l_seq=_seq.GetLength();
  unsigned l_target;

  unsigned start_index;
  if(_V_end+1<_flank_length)
    {//plus one because, we started from one past where the V aligned ends;
      start_index=0;
    }
  else
    {
      start_index=_V_end+1-_flank_length;
    }
    
  unsigned end_index=_J_start-1+_flank_length;
  if(end_index>=l_seq)
    {
      end_index=l_seq-1;
    }
  //cout<<"_J_start:"<<_J_start<<";flank_length:"<<_flank_length
  //    <<";end_index:"<<end_index<<endl;
  
  string seqNDN=_seq.GetSequence().substr(start_index, end_index-start_index+1);

  int flank_offset=start_index-1;//this one is one less than the start_index.

  _D.D_max_errors=(int)(AlignmentSettings::max_length_D_genes*_sm->GetScore('A','A')/(_sm->GetScore('A','A')-_sm->GetScore('A','C')));  
  //cout<<"-->D_max_errors:"<<_D.D_max_errors<<endl;
  //start doing the alignment by calling the localAlign function
  if(seqNDN.size()==0)
    {
      std::memset(_D.numOfAligned, 0, _numOfDSegs*sizeof(unsigned)/sizeof(char));
      return false;
    }
  
  //good, let's do the alignment
  SequenceString ss_seqNDN("temp", seqNDN);
  unsigned curr_numOfAligned;
  //cout<<"\t===>before initialize "<<endl;
  _D.initialize(_numOfDSegs);
  AlignmentString* las;
  SequenceString ss_target;
  //cout<<"\t====>bfore running the main for loop"<<endl;
  for(unsigned d=0;d<_numOfDSegs;d++)
    {
      //cout<<"\t****looping d:"<<d<<endl;
      ss_target=_genDs[d].Get_Seq();
      l_target=ss_target.GetLength();
      //cout<<"\t===>before calling local align"<<endl;
      //cout<<ss_seqNDN.toString()<<endl;
      //cout<<ss_target.toString()<<endl;
      LocalAlignment la(&ss_seqNDN, &ss_target, _sm,
			-100.0/*gapopen*/ ,-100.0 /*gapextension*/,
			1.0/*scale*/, _max_aligns, 0/*this is the gap model, 0 for affine model*/);
      //first determine the number of aligned for this current one
      //cout<<"\t****done with local alignment"<<endl;
      curr_numOfAligned=la.GetNumberOfAlignments();
      las=la.GetAlignmentArr();

      _D.numOfAligned[d]=curr_numOfAligned;
      _D.align_length[d]=new unsigned[curr_numOfAligned];
      _D.score[d]=new double[curr_numOfAligned];
      _D.n_errors[d]=new unsigned[curr_numOfAligned];
      _D.error_positions[d]=new unsigned*[curr_numOfAligned];
      unsigned* temp_error_positions;
      
      _D.align_position_left[d]=new unsigned [curr_numOfAligned];
      _D.align_position_right[d]=new unsigned [curr_numOfAligned];
      
      _D.deletions_left[d]=new unsigned [curr_numOfAligned];
      _D.deletions_right[d]=new unsigned [curr_numOfAligned];
      //cout<<"\t%%%%%before parsing the parms, aligned:"<<curr_numOfAligned<<endl;
      for(unsigned i=0;i<curr_numOfAligned;i++)
	{
	  //cout<<"\t\t@@@@@@@--->looping through aligned i:"<<i<<endl;
	  //cout<<"alignment #"<<i<<":"<<endl;
	  //cout<<las[i].toString();
	  _D.align_length[d][i]=las[i].GetPatternIndexEnd()-las[i].GetPatternIndexStart()+1;
	  _D.score[d][i]=las[i].GetScore();
	  if(flank_offset>=0)
	    {

	      _D.align_position_left[d][i]=las[i].GetPatternIndexStart()+flank_offset+1; //plus one, because flank_offset+1 is where the sequence started.
	      _D.align_position_right[d][i]=las[i].GetPatternIndexEnd()+flank_offset+1;
	    }
	  else
	    {
	      _D.align_position_left[d][i]=las[i].GetPatternIndexStart();
	      _D.align_position_right[d][i]=las[i].GetPatternIndexEnd();
	    }
	  _D.deletions_left[d][i]=las[i].GetSubjectIndexStart();
	  _D.deletions_right[d][i]=_genDs[d].Get_Seq().GetLength()-las[i].GetSubjectIndexEnd()-1;

	  /*cout<<"\t\talign_length:"<<_D.align_length[d][i]<<endl;
	  cout<<"\t\tscore:"<<_D.score[d][i]<<endl;
	  cout<<"\t\talign_position_left:"<<_D.align_position_left[d][i]<<endl;
	  cout<<"\t\talign_position_right:"<<_D.align_position_right[d][i]<<endl;
	  cout<<"\t\tdeletions_left:"<<_D.deletions_left[d][i]<<endl;
	  cout<<"\t\tdeletions_right:"<<_D.deletions_right[d][i]<<endl;
	  */

	  temp_error_positions=new unsigned[_D.D_max_errors];
	  //cout<<"\t\t****beofre finding error fucntions"<<endl;
	  unsigned pos_start1=_D.align_position_left[d][i]-flank_offset-1;
	  unsigned pos_end1=_D.align_position_right[d][i]-flank_offset-1;
	  if(flank_offset<0)
	    {
	      pos_start1=_D.align_position_left[d][i];
	      pos_end1=_D.align_position_right[d][i];
	    }
	  _D.n_errors[d][i]=findErrors
	    (seqNDN, _genDs[d].Get_Sequence(), 
	     pos_start1, pos_end1,
	     _D.deletions_left[d][i],l_target-_D.deletions_right[d][i]-1,
	     _D.D_max_errors,flank_offset, temp_error_positions);

	  //cout<<"\t\tafter finding the errors"<<endl;
	  //now copy over the elements
	  _D.error_positions[d][i]=new unsigned[_D.n_errors[d][i]];
	  memcpy(_D.error_positions[d][i], temp_error_positions,
		 _D.n_errors[d][i]*sizeof(unsigned)/sizeof(char));
	  //cout<<"\t\tcleaning up mem"<<endl;
	  //clean up the memory
	  delete [] temp_error_positions;
	  //cout<<"\t\tdone!!!"<<endl;
	}
    }//end of num of D gene segs for loop
  
  //now we are ready with alignment, first sort in order to figure
  //out the D.allele_order, prepare the index array first
  //cout<<"\n===>sort to get the allele_order "<<endl;
  double* highestScore=new double[_numOfDSegs];
  //cout<<"before sorting score array:"<<endl;
  for(unsigned i=0;i<_numOfDSegs;i++)
    {
      _D.allele_order[i]=i;
      highestScore[i]=_D.score[i][0];
      //cout<<highestScore[i]<<",";
    }
  //cout<<endl;
  QuickSort(highestScore, 0, _numOfDSegs-1, _D.allele_order,NULL);
  //cout<<"after sorting score array:"<<endl;
  /*for(unsigned i=0;i<_numOfDSegs;i++)
    {
      //_D.allele_order[i]=i;
      //highestScore[i]=_D.score[i][0];
      cout<<highestScore[i]<<",";
    }
  cout<<endl;
  */
  Reverse(_D.allele_order, _numOfDSegs);

  //cout<<"after sorting index array:"<<endl;
  /*for(unsigned i=0;i<_numOfDSegs;i++)
    {
      //_D.allele_order[i]=i;
      //highestScore[i]=_D.score[i][0];
      cout<<_D.allele_order[i]<<",";
    }
  cout<<endl;
  */
  //cout<<"\t****done with sorting"<<endl;

  //now we are ready to take care of p_nucleotides and negative excess error
  //first need to initialize the p_region array to all zeros
  //cout<<"\n----->>>>>Ready to figure out p_region and excess error"<<endl;
  _D.p_region_max_length_left=new unsigned**[_numOfDSegs];
  _D.p_region_max_length_right=new unsigned**[_numOfDSegs];
  _D.excess_error_positions_left=new unsigned**[_numOfDSegs];
  _D.excess_error_positions_right=new unsigned**[_numOfDSegs];

  //cout<<"\tbefore looping......_D_maximum_deletion:"<<_D_maximum_deletion<<endl;
  for(unsigned i=0;i<_numOfDSegs;i++)
    {
      //cout<<"\t====>looping....i:"<<i<<endl;
      _D.p_region_max_length_left[i]=new unsigned*[_D.numOfAligned[i]];
      _D.p_region_max_length_right[i]=new unsigned*[_D.numOfAligned[i]];
      _D.excess_error_positions_left[i]=new unsigned*[_D.numOfAligned[i]];
      _D.excess_error_positions_right[i]=new unsigned*[_D.numOfAligned[i]];

      for(unsigned j=0;j<_D.numOfAligned[i];j++)
	{
	  _D.p_region_max_length_left[i][j]=
	    new unsigned[1+_D_maximum_deletion];
	  std::memset(_D.p_region_max_length_left[i][j],0,
		      (1+_D_maximum_deletion)*sizeof(unsigned)/sizeof(char));
	  _D.p_region_max_length_right[i][j]=
	    new unsigned[1+_D_maximum_deletion];
	  std::memset(_D.p_region_max_length_right[i][j],0,
		      (1+_D_maximum_deletion)*sizeof(unsigned)/sizeof(char));
	  _D.excess_error_positions_left[i][j]=
	    new unsigned[_negative_excess_deletion_max];
	  std::memset(_D.excess_error_positions_left[i][j],0,
		      _negative_excess_deletion_max*sizeof(unsigned)/sizeof(char));
	  _D.excess_error_positions_right[i][j]=
	    new unsigned[_negative_excess_deletion_max];
	  std::memset(_D.excess_error_positions_right[i][j],0,
		      _negative_excess_deletion_max*sizeof(unsigned)/sizeof(char));
	}
    }
  //cout<<"\t\t*****done initializing the arrays"<<endl;
  //now ready to call to find p-nucleotide and excess error
  DeterminePalindromAndExcessError_D
    (_seq, _genDs, _D.deletions_left, _D.deletions_right, 
     _negative_excess_deletion_max, _D_maximum_deletion,
     _D.align_length, _numOfDSegs, _D.numOfAligned,
     _D.align_position_left, _D.align_position_right,
     _D.p_region_max_length_left, _D.p_region_max_length_right,
     _D.excess_error_positions_left, 
     _D.excess_error_positions_right
     );
  //cout<<"@@@@@@@@finally done!!!!"<<endl;
  //done!!!
  return true;
}//end of function MatchD

/*find the error positions for two sequence aligned.
 * the output array has to be allocated by the caller outside 
 * NOTE: BE CAREFUL here, the output array error_positions
 * holding the absolute position/index, starting at the 
 * beginning of the sequence _seq1. NOT a relative position
 * flank_offset is one less where the seqNDN starts, seqNDN is the subSting
 *        we want to compare to do the alignment. we have this flank_offset
 *        in order to get the absoluted error_position instead a relative one
 */
unsigned findErrors(const string& _seq1, const string& _seq2,
	   const unsigned& pos_start1, const unsigned& pos_end1,
		    const unsigned& pos_start2, const unsigned& pos_end2,
		    const unsigned& max_n_errors, const int& _flank_offset,
		    /*output*/unsigned* error_position
		    )
{
  //first, do some checking
  if(pos_end1-pos_start1!=pos_end2-pos_start2)
    {
      throw "bad array index, the seq to be check for match/mismatch are not of same length!";
    }
  
  unsigned n_error=0;
  for(unsigned i=0;i<pos_end2-pos_start2+1;i++)
    {
      if(_seq1.at(pos_start1+i)!=_seq2.at(pos_start2+i))
	{	  
	  if(n_error<max_n_errors)
	    {
	      if(_flank_offset>=0)
		error_position[n_error]=i+pos_start1+_flank_offset+1;
	      else
		error_position[n_error]=i+pos_start1;
	    }
	  n_error++;
	}
    }
  
  return n_error;
}

/*Again in this function, the output has been correctly initialized
 *correctly in the caller outside.
 */
void DeterminePalindromAndExcessError_D
( const SequenceString& _seq, const GenomicD* _genDs,
  /*const unsigned* _ok_order,*/ unsigned** _deletions_left,
  unsigned** _deletions_right,
  const unsigned& _negative_excess_deletions_max, 
  const unsigned& _D_maximum_deletion,
  const unsigned* const* _align_length, const unsigned& _numOfDSegs,
  const unsigned* _numOfAligned, 
  const unsigned* const* _align_position_left, const unsigned* const* _align_position_right,
  /*output*/ unsigned*** _p_region_max_length_left, 
  unsigned*** _p_region_max_length_right,
  unsigned*** _excess_error_positions_left,
  unsigned*** _excess_error_positions_right
  )
{
  //cout<<"++++>inside determine palindrom function"<<endl;
  bool still_palindrome=true;
  unsigned p=0;
  unsigned l_seq=_seq.GetLength();
  unsigned l_target;
  //now go through to find palindromic nucleotides for each alignment
  for(unsigned d=0;d<_numOfDSegs;d++)
    {
      //cout<<"\tmain loop:"<<d<<endl;
      string target=_genDs[d].Get_Sequence();
      l_target=target.size();
      //for each aligned in this D seg
      for(unsigned na=0;na<_numOfAligned[d];na++)
	{
	  //cout<<"\t\tsubmain loop:"<<na<<"/"<<_numOfAligned[d]<<endl;
	  //getting the left p-nucleotides first
	  //for each possible deletions
	  int nd=_deletions_left[d][na]-_negative_excess_deletions_max;
	  if(nd<0)
	    nd=0;
	  int max_nd=_D_maximum_deletion;
	  if(max_nd>(signed)(_align_length[d][na]+_deletions_left[d][na]))
	    max_nd=_align_length[d][na]+_deletions_left[d][na];
	  //cout<<"\t\t==???max_nd:"<<max_nd<<endl;
	  for(;nd<=max_nd;nd++)
	    {
	      //cout<<"\t\t\tloop nd:"<<nd<<endl;
	      //for each value of deletions, find longest half-p from the implied end of the gene sequence
	      p=0;
	      still_palindrome=true;
	      //upper bound on p is the length of the match seq or seq to the left of the aligned ones
	      int max_p_length=_align_length[d][na]-(nd-_deletions_left[d][na]);
	      if(max_p_length<0)
		max_p_length=0;
	      int tempV=_align_position_left[d][na]+(nd-_deletions_left[d][na]);
	      if(max_p_length>tempV)
		{
		  max_p_length=tempV;
		}
	      //cout<<"max_p_length :"<<max_p_length<<endl;
	      while(still_palindrome&&((signed)p<(signed)max_p_length))
		{
		  still_palindrome=target.at(nd+p)==DnaComplement(_seq.GetSequence().at(_align_position_left[d][na]-p+nd-_deletions_left[d][na]-1));
		  if(still_palindrome)
		    {
		      p++;
		    }
		}
	      //cout<<"\t\t\tset the value nd:"<<nd<<endl;
	      _p_region_max_length_left[d][na][nd]=p;
	    }//end of for nd<max_nd
	  //cout<<"\t\t==*****>end of left palindrom"<<endl;
	  //start doing excess error 
	  unsigned tempArray[]={_negative_excess_deletions_max, _deletions_left[d][na], _align_position_left[d][na]};
	  unsigned n_excess=min_mf(tempArray,3);
	  unsigned k;
	  unsigned runningIndex_excessError=0;
	  //find the mismatchs
	  for(k=0;k<n_excess;k++)
	    {
	      if(target.at(_deletions_left[d][na]-n_excess+k)!=
		 _seq.GetSequence().at(_align_position_left[d][na]-n_excess+k))
		{
		  _excess_error_positions_left[d][na][runningIndex_excessError]=
		    _align_position_left[d][na]-n_excess+k;
		  runningIndex_excessError++;
		}
	    }//end n_excess
	  for(;runningIndex_excessError<_negative_excess_deletions_max;runningIndex_excessError++)
	    {
	      _excess_error_positions_left[d][na][runningIndex_excessError]=0;
	    }
	  //cout<<"\t\tend of left excess error"<<endl;
	  
	  //===>>>>doing right side things
	  nd=_deletions_right[d][na]-_negative_excess_deletions_max;
	  if(nd<0)
	    nd=0;
	  max_nd=_D_maximum_deletion;
	  if(max_nd>(signed)(_align_length[d][na]+_deletions_right[d][na]))
	    max_nd=_align_length[d][na]+_deletions_right[d][na];
	  	  
	  for(;nd<=max_nd;nd++)
	    {
	      //cout<<"\t\tloop nd:"<<nd<<endl;
	      //for each value of deletions, find longest half-p from the implied end of the gene sequence
	      p=0;
	      still_palindrome=true;
	      //upper bound on p is the length of the match seq or seq to the left of the aligned ones
	      int max_p_length=_align_length[d][na]-(nd-_deletions_right[d][na]);
	      if(max_p_length<0)
		max_p_length=0;
	      //cout<<"right bound:"<<l_seq-_align_position_right[d][na]+(nd-_deletions_right[d][na])<<endl;
	      int tempV=l_seq-1-_align_position_right[d][na]+(nd-_deletions_right[d][na]);
	      if(max_p_length>tempV)
		{
		  max_p_length=tempV;
		}

	      while(still_palindrome&&((signed)p<(signed)max_p_length))
		{
		  //cout<<"while loop "<<p<<endl;
		  still_palindrome=target.at(l_target-nd-p-1)==DnaComplement(_seq.GetSequence().at(_align_position_right[d][na]+1+p-(nd-_deletions_right[d][na])));
		  if(still_palindrome)
		    {
		      p++;
		    }
		}
	      //cout<<"set p value:"<<p<<endl;
	      _p_region_max_length_right[d][na][nd]=p;
	    }//end of for nd<max_nd
	  //cout<<"\t\tend of right palindrom"<<endl;
	  
	  //start doing the excess error for right side.
	  //tempArray[0]=_negative_excess_deletions_max;
	  tempArray[1]=_deletions_right[d][na]; 
	  tempArray[2]=l_seq-_align_position_right[d][na]-1;
      
	  n_excess=min_mf(tempArray,3);
	  //unsigned k;
	  
	  runningIndex_excessError=0;
	  //find the mismatchs
	  for(k=0;k<n_excess;k++)
	    {
	      if(target.at(l_target-_deletions_right[d][na]+k)!=
		 _seq.GetSequence().at(_align_position_right[d][na]+1+k))
		{
		  _excess_error_positions_right[d][na][runningIndex_excessError]=
		    _align_position_right[d][na]+1+k;
		  runningIndex_excessError++;
		}
	    }//end n_excess
	  for(;runningIndex_excessError<_negative_excess_deletions_max;runningIndex_excessError++)
	    {
	      _excess_error_positions_right[d][na][runningIndex_excessError]=0;
	    }
	  //done!!!
	  //cout<<"\t\tend of right excess error"<<endl;

	}//end of _numOfAligned
    }//end of outer for loop d : numOfDSegs
}//end of function find p-nucleotides

//usually, D alignment is successful anyway.
void Alignment_D::Serialize(ofstream& _ofs)const
{
  //do necessary checking
  if(!_ofs.is_open())
    {
      cout<<"**ERROR**: closed file buffer. quit..."<<endl;
      exit(-1);
    }
  //serializing......
  char* p_char; //pointer used to direct the writing to the file stream
  //save n_D_alleles 
  p_char=(char*)&n_D_alleles;
  _ofs.write(p_char, sizeof(unsigned));
  
  //D_max_errors
  p_char=(char*)&D_max_errors;
  _ofs.write(p_char, sizeof(unsigned));
  
  //numOfAligned --- in here, we don't check for successful alignment, 
  //because align d is usually sucessful.
  p_char=(char*)numOfAligned;
  _ofs.write(p_char, sizeof(unsigned)*n_D_alleles);

  //align_length  
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)align_length[i];
      _ofs.write(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //score
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)score[i];
      _ofs.write(p_char, sizeof(double)*numOfAligned[i]);
    }

  //n_errors
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)n_errors[i];
      _ofs.write(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //error_positions
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      for(unsigned j=0;j<numOfAligned[i];j++)
	{	  
	  p_char=(char*)error_positions[i][j];
	  _ofs.write(p_char, sizeof(unsigned)*n_errors[i][j]);
	}
    }
  
  //align_position_left
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)align_position_left[i];
      _ofs.write(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //align_position_right
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)align_position_right[i];
      _ofs.write(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //deletion_left
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)deletions_left[i];
      _ofs.write(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //deletion_right
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_char=(char*)deletions_right[i];
      _ofs.write(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //p_region_max_length_left
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      for(unsigned j=0;j<numOfAligned[i];j++)
	{	  
	  p_char=(char*)p_region_max_length_left[i][j];
	  _ofs.write(p_char, sizeof(unsigned)*(1+AlignmentSettings::D_maximum_deletion));
	}
    }
  
  //p_region_max_length_right
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      for(unsigned j=0;j<numOfAligned[i];j++)
	{	  
	  p_char=(char*)p_region_max_length_right[i][j];
	  _ofs.write(p_char, sizeof(unsigned)*(1+AlignmentSettings::D_maximum_deletion));
	}
    }
  
  //excess_error_position
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      for(unsigned j=0;j<numOfAligned[i];j++)
	{	  
	  p_char=(char*)excess_error_positions_left[i][j];
	  _ofs.write(p_char, sizeof(unsigned)*(AlignmentSettings::negative_excess_deletions_max));
	}
    }
  
  //excess_error_position
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      for(unsigned j=0;j<numOfAligned[i];j++)
	{	  
	  p_char=(char*)excess_error_positions_right[i][j];
	  _ofs.write(p_char, sizeof(unsigned)*(AlignmentSettings::negative_excess_deletions_max));
	}
    }
    
  //allele_order
  p_char=(char*)allele_order;
  _ofs.write(p_char, sizeof(unsigned)*n_D_alleles);

  //done  
}

//usually, D alignment is anyway sucessfully
//NOTE: here in side this function, we assume all fields are the default,
//zero or NULL.
void Alignment_D::Deserialize(ifstream& _ifs)
{
  //cout<<"inside the deserialize function"<<endl;
  //first we need to check for there is fields to read
  if(_ifs.eof())
    {
      cout<<"*****ERROR: can not read fields, end of the file"<<endl;
      exit(-1);
    }

  //now do the reading.
  // one important thing here is we want to check whether the original
  //object is empty, it could be nonempty, so we want to clean it
  //up while we are reading the new one
  //cout<<"n_D_alleles:"<<n_D_alleles<<endl;
  unsigned original_n_D_alleles=n_D_alleles;
  //unsigned original_max_errors=D_max_errors;
  unsigned * original_numOfAligned;
  if((signed)original_n_D_alleles!=0)
    original_numOfAligned=new unsigned[original_n_D_alleles];
  else
    original_numOfAligned=new unsigned[34];//AlignmentSettings::n_D_alleles];
  //read n_D_alleles
  char* p_char=(char*)&n_D_alleles;
  _ifs.read(p_char, sizeof(unsigned));

  //max_error
  p_char =(char*)&D_max_errors;
  _ifs.read(p_char, sizeof(unsigned));
  
  //numOfaligned
  if(numOfAligned!=NULL)
    {
      //cout<<"\tcopy over the non-null numOfAligned"<<endl;
      memcpy(original_numOfAligned, numOfAligned, sizeof(unsigned)*original_n_D_alleles);
      delete[] align_length;
    }

  numOfAligned=new unsigned [n_D_alleles];
  p_char=(char*)numOfAligned;
  _ifs.read(p_char, sizeof(unsigned)*n_D_alleles);

  //align_length
  //cout<<"align_length"<<endl;
  if(align_length!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(align_length[i]!=NULL)
	    delete [] align_length[i];
	}
      delete [] align_length;
    }
  align_length=new unsigned* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      align_length[i]=new unsigned[numOfAligned[i]];
      p_char=(char*)align_length[i];
      _ifs.read(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //score
  //cout<<"score"<<endl;
  if(score!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(score[i]!=NULL)
	    delete [] score[i];
	}
      delete [] score;
    }
  score=new double* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      score[i]=new double[numOfAligned[i]];
      p_char=(char*)score[i];
      _ifs.read(p_char, sizeof(double)*numOfAligned[i]);
    }

  //n_errors
  if(n_errors!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(n_errors[i]!=NULL)
	    delete [] n_errors[i];
	}
      delete [] n_errors;
    }
  n_errors=new unsigned* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      n_errors[i]=new unsigned[numOfAligned[i]];
      p_char=(char*)n_errors[i];
      _ifs.read(p_char, sizeof(unsigned)*numOfAligned[i]);
    }

  //error_positions
  //cout<<"error_positions"<<endl;
  if(error_positions!=NULL)
    {
      //cout<<"non-null value"<<endl;
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(error_positions[i]!=NULL)
	    {
	      for(unsigned j=0;j<original_numOfAligned[i];j++)
		{
		  if(error_positions[i][j]!=NULL)
		    {
		      delete [] error_positions[i][j];
		    }
		}
	      delete [] error_positions[i];
	    }
	}
      delete [] error_positions;
    }
  //cout<<"alloc **pointer"<<endl;
  error_positions=new unsigned** [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      //cout<<"\t==>alloc *pointer, numOfAligned[i]:"<<numOfAligned[i]<<endl;
      error_positions[i]=new unsigned*[numOfAligned[i]];
      for(unsigned j=0;j<numOfAligned[i];j++)
	{
	  //cout<<"\t\talloc pointer, n_errors[i][j]:"<<n_errors[i][j]<<endl;
	  error_positions[i][j]=new unsigned[n_errors[i][j]];
	  p_char=(char*)error_positions[i][j];
	  _ifs.read(p_char, sizeof(unsigned)*n_errors[i][j]);
	}
    }

  //align_position_left
  //cout<<"aign_position left"<<endl;
  if(align_position_left!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(align_position_left[i]!=NULL)
	    delete [] align_position_left[i];
	}
      delete [] align_position_left;
    }
  align_position_left=new unsigned* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      align_position_left[i]=new unsigned[numOfAligned[i]];
      p_char=(char*)align_position_left[i];
      _ifs.read(p_char, sizeof(unsigned)*numOfAligned[i]);
    }
  //align_position_right
  if(align_position_right!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(align_position_right[i]!=NULL)
	    delete [] align_position_right[i];
	}
      delete [] align_position_right;
    }
  align_position_right=new unsigned* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      align_position_right[i]=new unsigned[numOfAligned[i]];
      p_char=(char*)align_position_right[i];
      _ifs.read(p_char, sizeof(unsigned)*numOfAligned[i]);
    }
  //deletion_left
  //cout<<"deletion left"<<endl;
  if(deletions_left!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(deletions_left[i]!=NULL)
	    delete [] deletions_left[i];
	}
      delete [] deletions_left;
    }
  deletions_left=new unsigned* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      deletions_left[i]=new unsigned[numOfAligned[i]];
      p_char=(char*)deletions_left[i];
      _ifs.read(p_char, sizeof(unsigned)*numOfAligned[i]);
    }
  //deletion_right
  if(deletions_right!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(deletions_right[i]!=NULL)
	    delete [] deletions_right[i];
	}
      delete [] deletions_right;
    }
  deletions_right=new unsigned* [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      deletions_right[i]=new unsigned[numOfAligned[i]];
      p_char=(char*)deletions_right[i];
      _ifs.read(p_char, sizeof(unsigned)*numOfAligned[i]);
    }
  
  //p_region_max_length_left
  if(p_region_max_length_left!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(p_region_max_length_left[i]!=NULL)
	    {
	      for(unsigned j=0;j<original_numOfAligned[i];j++)
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
    }
  
  p_region_max_length_left=new unsigned** [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_region_max_length_left[i]=new unsigned*[numOfAligned[i]];
      for(unsigned j=0;j<numOfAligned[i];j++)
	{
	  p_region_max_length_left[i][j]
	    =new unsigned[AlignmentSettings::D_maximum_deletion+1];
	  p_char=(char*)p_region_max_length_left[i][j];
	  _ifs.read(p_char, sizeof(unsigned)*(1+AlignmentSettings::D_maximum_deletion));
	}
    }

  //p_region_max_length_right
  if(p_region_max_length_right!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(p_region_max_length_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<original_numOfAligned[i];j++)
		{
		  if(p_region_max_length_right[i][j]!=NULL)
		    {
		      delete [] p_region_max_length_right[i][j];
		    }
		}
	      delete [] p_region_max_length_right[i];
	    }
	}
      delete [] p_region_max_length_right;
    }
  
  p_region_max_length_right=new unsigned** [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      p_region_max_length_right[i]=new unsigned*[numOfAligned[i]];
      for(unsigned j=0;j<numOfAligned[i];j++)
	{
	  p_region_max_length_right[i][j]
	    =new unsigned[AlignmentSettings::D_maximum_deletion+1];
	  p_char=(char*)p_region_max_length_right[i][j];
	  _ifs.read(p_char, sizeof(unsigned)*(1+AlignmentSettings::D_maximum_deletion));
	}
    }
  //excess_error_position_left
  if(excess_error_positions_left!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(excess_error_positions_left[i]!=NULL)
	    {
	      for(unsigned j=0;j<original_numOfAligned[i];j++)
		{
		  if(excess_error_positions_left[i][j]!=NULL)
		    {
		      delete [] excess_error_positions_left[i][j];
		    }
		}
	      delete [] excess_error_positions_left[i];
	    }
	}
      delete [] excess_error_positions_left;
    }
  
  excess_error_positions_left=new unsigned** [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      excess_error_positions_left[i]=new unsigned*[numOfAligned[i]];
      for(unsigned j=0;j<numOfAligned[i];j++)
	{
	  excess_error_positions_left[i][j]
	    =new unsigned[AlignmentSettings::negative_excess_deletions_max];
	  p_char=(char*)excess_error_positions_left[i][j];
	  _ifs.read(p_char, sizeof(unsigned)*(AlignmentSettings::negative_excess_deletions_max));
	}
    }  
  
  //excess_error_position_right
  if(excess_error_positions_right!=NULL)
    {
      for(unsigned i=0;i<original_n_D_alleles;i++)
	{
	  if(excess_error_positions_right[i]!=NULL)
	    {
	      for(unsigned j=0;j<original_numOfAligned[i];j++)
		{
		  if(excess_error_positions_right[i][j]!=NULL)
		    {
		      delete [] excess_error_positions_right[i][j];
		    }
		}
	      delete [] excess_error_positions_right[i];
	    }
	}
      delete [] excess_error_positions_right;
    }
  
  excess_error_positions_right=new unsigned** [n_D_alleles];
  for(unsigned i=0;i<n_D_alleles;i++)
    {
      excess_error_positions_right[i]=new unsigned*[numOfAligned[i]];
      for(unsigned j=0;j<numOfAligned[i];j++)
	{
	  excess_error_positions_right[i][j]
	    =new unsigned[AlignmentSettings::negative_excess_deletions_max];
	  p_char=(char*)excess_error_positions_right[i][j];
	  _ifs.read(p_char, sizeof(unsigned)*(AlignmentSettings::negative_excess_deletions_max));
	}
    }
  //allele_order
  if(allele_order!=NULL)
    {
      delete [] allele_order;
    }
  allele_order=new unsigned [n_D_alleles];
  p_char=(char*)allele_order;
  _ifs.read(p_char, sizeof(unsigned)*n_D_alleles);
    
  //clean up memory
  delete [] original_numOfAligned;
}


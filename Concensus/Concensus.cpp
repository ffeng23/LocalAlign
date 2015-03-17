#include "Concensus.hpp"
#include "../SIGPIG/MatrixFunctions.hpp"

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

//run the concensus for each individual file, which will give out one concensus sequence
//
SequenceString GetConcensus(const string& _fileName, unsigned& _numOfSeq)
{
  srand(time(NULL));

  cout<<"---->doing file:"<<_fileName<<endl;
  //calling the FastaHandler to handle it
  vector<SequenceString> seq;
  unsigned totalNumber=ReadFasta(_fileName, seq, true);
  _numOfSeq=totalNumber;
  SequenceString c_seq; //consensus seq
  c_seq.SetName(_fileName.substr(0, _fileName.length()-5));
  
  if(totalNumber<=0)
    return c_seq;
  //now go through each position to figure out the concensus
  //first get the length (max)of the string
  unsigned len=seq.at(0).GetLength();
  for(unsigned i=i;i<totalNumber;i++)
    {
      if(len<seq.at(i).GetLength())
	len=seq.at(i).GetLength();
    }
  //cout<<"len: "<<len<<endl;
  //now we got the length of the consensus
  //we assume the sequences in the input file are aligned well
  //they could be of different length, but should be fixed in the front
  
  //c_seq.SetSequence(string(len, 'A'));
  char nt[4]={'A', 'C','G', 'T'};
  char currentLetter;
  unsigned stats[5]={0,0,0,0,0}; 
  unsigned int stats_index[5]={0,1,2,3,4};
  //remembering the count stats for each "A/a, C/c, G/g, T/t or something else"
  
  //now go through to get the concensus

  string sequence_con(len, 'A'); //set the concensus sequence to be all 'A' in the beginning.
  //cout<<"start doing the loop"<<endl;
  for(unsigned i=0;i<len;i++)
    {
      //cout<<"loop i:"<<i<<endl;
      stats[0]=0;stats[1]=0;stats[2]=0;stats[3]=0;stats[4]=0;
      stats_index[0]=0;stats_index[1]=1;stats_index[2]=2;stats_index[3]=3;stats_index[4]=4;
      for(unsigned j=0;j<totalNumber;j++)
	{
	  if(seq.at(j).GetLength()<=i) //this current one is too short, we go next
	    continue;
	  switch(seq.at(j).GetSequence().at(i))
	    {
	    case 'A':
	    case 'a':
	      stats[0]++;
	      break;
	    case 'C':
	    case 'c':
	      stats[1]++;
	      break;
	    case 'G':
	    case 'g':
	      stats[2]++;
	      break;
	    case 'T':
	    case 't':
	      stats[3]++;
	      break;
	    default:
	      stats[4]++;
	      cout<<"WARNING:found a nt that is undetermined (not A/C/G/T)"<<endl;
	      cerr<<"WARNING:found a nt that is undetermined (not A/C/G/T)"<<endl;
	      break;
	    
	    }//end of switch
	  
	}//end of inner for all sequences
      //now for "this" position, we need to get the concensus
      //check for the best stats
      //cout<<"do sorting...."<<endl;
      QuickSort<unsigned>(stats, 0,5,stats_index,NULL);//in ascending order
      //stats now is in ascending order

      unsigned bestStats
	=stats[4];//biggest one
      //check for ties
      unsigned numberOfTies=0;
      for(unsigned k=0;k<4;k++)
	{
	  if(bestStats==stats[k])
	    {
	      numberOfTies++;
	    }
	}
      //cout<<"numberOfTies:"<<numberOfTies<<endl;
      
      //if numberOfTies is not zero, means we have ties,
      //need to decide which one to get
      //run random number to get randomly picked one
      if(numberOfTies==0)
	{
	  //no tie, pick the best one
	  if(stats_index[4]!=4) //this is a good one
	    {
	      currentLetter=nt[stats_index[4]];
	    }
	  else //the best one is unknow letter, we need to pick the next best
	    {
	      cout<<"WARNING: best one is the unknow letter, try random select"<<endl;
	      currentLetter=nt[stats_index[3]];
	    }

	}
      else //we get a tie,
	{
	  //randomly pick one from the ties
	  //cout<<"===tie breaker now"<<endl;
	  unsigned index_r=4;
	  while(index_r==4)
	    {
	      //cout<<"keep working...."<<endl;
	      index_r=stats_index[4-rand()%(numberOfTies+1)];
	    }
	  currentLetter=nt[index_r];
	}
      //now set the letter to the c_seq
      sequence_con[i]=currentLetter;
      /*if(i%50==0)
	{
	  cout<<".";
	  fflush(stdout);
	  }*/
    }//end of positions of c_seq
  //now we are doine.
  cout<<"done!!"<<endl;
  c_seq.SetSequence(sequence_con);
  return c_seq;
}


//go through a bunch of input file and generate one sequence from each file,
//we also want to remember how many sequences for each concensus sequence.
void GenerateConcensus(const string* _ifileNames, const unsigned& _numOfFiles,
		       /*output*/ SequenceString* _ss, unsigned* countOfConcensus
		       )
{
  //for each file, run the concensus 
  for(unsigned i=0;i<_numOfFiles;i++)
    {
      _ss[i]=GetConcensus(_ifileNames[i], countOfConcensus[i]);
    }
 }

//based one the concensus sequences, write up the input sequence file for our analys
//it has specific format, like 
//header
//sequence count
//see the ../sample.data for the detail format.
//this format is adopted from the Matlab code.
void GenerateSequenceFile(const SequenceString* _ss, const unsigned* countOfConcensus,
			  const unsigned& _numOfSeq, const string& _ofname)
  {
    ofstream ofs_p(_ofname.c_str());
      
  if(!ofs_p.is_open())
    {
      cout<<">>>>>>ERROR:the input file \""<<_ofname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  //start writing the header
  ofs_p<<"#Pool=mem\n"
       <<"#Frame=coding\n"
       <<"#N_Sequences="<<_numOfSeq<<"\n"
       <<"#Read_Length=-1\n"
       <<"#N_Reads="<<_numOfSeq<<"\n"
       <<"#N_Reads_Pool="<<_numOfSeq<<"\n"
       <<"Sequence\tCounts\n";
  for(unsigned i=0;i<_numOfSeq;i++)
    {
      ofs_p<<_ss[i].GetSequence()<<"\t"<<countOfConcensus[i]<<"\n";
    }
  //done
  ofs_p.close();
  }


//this is the function to put together everything, so the main function only
//providing the directory and this function will go through it and generate concensus
//write a sequence input file.
//
bool DoGenerateSequenceFile(const char* _path, const string& _ofname)
{
  //calling first to get the file array
  string* files=NULL;
  unsigned numOfFiles=0;
  cout<<"Collecting the files with \".fasta\" formats...."<<endl;
  GetFileNames(_path, &files, numOfFiles);
  cout<<"Done....."<<numOfFiles<<" found!"<<endl;

  cout<<"Start generating concensus...."<<endl;
  SequenceString* ss=new SequenceString[numOfFiles];
  unsigned* countOfConcensus=new unsigned[numOfFiles];
  GenerateConcensus(files, numOfFiles, ss, countOfConcensus);
  cout<<"Done with concensus!"<<endl;
  //
  cout<<"Start writing the seuqnece files....."<<endl;
  GenerateSequenceFile(ss, countOfConcensus, numOfFiles, _ofname);
  cout<<"Done!!"<<endl;
  //mem
  delete[] files;
  delete[] ss;
  delete[] countOfConcensus;
		    
  return true;
}

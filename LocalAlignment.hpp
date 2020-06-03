#ifndef LOCALALIGNMENT_HPP
#define LOCALALIGNMENT_HPP
#include "pairwiseAlignment.hpp"
#include <vector>

using namespace std;

//this entry class is used to define each local alignment entry/path
class Path
{
public:
  Path();  //~Path() in this case we don't need a destructor since we don't have dynamic allocation of memeory.
  Path( const unsigned int _OptimalIndex[2], const unsigned int _startIndex [2], const double& _optimalValue );
  void SetOptimalIndex(unsigned int _OptimalIndex[2]);
  void SetStartIndex(unsigned int _startIndex [2]);
  void SetOptimalValue(const double& _optimalValue);
  const unsigned int* GetOptimalIndex() const ;
  const unsigned int* GetStartIndex() const;
  double GetOptimalValue() const;
  bool isTraced() const;
  void setTraced() ;
private:
  unsigned int c_optimalIndex[2];
  unsigned int c_startIndex[2];
  double c_optimalValue;
  bool c_traced;
};

class PathElementEntry
{
	public:
	PathElementEntry();
	~PathElementEntry();
	
	//copy constructor
	PathElementEntry(const PathElementEntry& pe);
	
	unsigned GetNumberOfPathes() const;
	unsigned GetPathID(unsigned index) const;
	double GetPathOptimalScore(unsigned index) const;
	unsigned GetPathOptimalIndexPattern(unsigned index) const;
	unsigned GetPathOptimalIndexSubject(unsigned index) const;
	//in here, we have to set up thing at the same time together
	// input index: which path this is for the entry in the pathElementEntry vector, like first, second, third....
	//			pathID: the path ID in the path vector
	//		optimalScore: the m
	void AddPathInfoNoCheckingDuplication(const unsigned& index, const unsigned& pathID, 
			const double& optimalScore, 
			const unsigned& optimalIndexPattern, const unsigned& optimalIndexSubject);
	//combining the extra one with this current one.
	void AddPathInfo(const PathElementEntry& pe);
	
	//remove a path (indicated by "index") from the records.
	void RemovePath(const unsigned index);
	//based on the pathid, look up the index of the path in the records.
	unsigned LookUpPathIndex(const unsigned& pathID) const ;

	string toString() const;
	private:
	PathElementEntry& operator = (const PathElementEntry& p) { return *this;}; //disable it  
	//pathID and c_optimalScore a pair to "rember" the information about
	//pathes the current entry engaged. each entry can be joining multiple
	//path, thanking to the cases where mismatch score and gap penalty are 
	//identical. it can be more one, more than 2 (might be rare)
	//unsigned c_pathNumber; //this is used to indicate how many pathes the current entry engaged.
	vector<unsigned> c_pathID;
	vector<double> c_optimalScore;//used to indicate up to this current entry the optimal score for the specific path
	vector<unsigned> c_optimalIndexPattern;//used to indicate where the optimal score on pattern is.
	vector<unsigned> c_optimalIndexSubject;//same thing, but on subject.
};


class LocalAlignment: public PairwiseAlignment
{
public:
  LocalAlignment(SequenceString* _pattern, SequenceString* _subject, 
		 const ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		 const double& _gextension=-5, const double& _scale=1,const int& _numOfAlignments=1, const short& _typeOfGapModel=1);//here we default to 1 454 markov chain model,
  
  virtual ~LocalAlignment();

  double* GetScoreArr();
  AlignmentString* GetAlignmentArr();
  unsigned int GetNumberOfAlignments();
  //void alignLM();
  
  void printPathEntryTable() const;
protected:
  LocalAlignment(); //default one is disabled, since we need to the information for pattern and subject length for
					//for deep destruction 
  virtual void align();
  virtual void traceBack();
	
	//in this function, we check the current node that is being used by the current 
	//path and then look for the other pathes that is going through this nodes.
	//here we simply write them down to a temporary vector for later update. By later
	//we mean after trace back this current path.
	//we set it up as a object function, since we need to get access to the information
	//about the class objec, the c_pattern length and the pathElementTable. 
	void recordPathInfoToBeRemove(const unsigned& _patternIndex, 
				const unsigned & _subjectIndex
			, const unsigned& c_current_path_tb,
					//outputs 
				vector<unsigned>& vec_pathToBeChecked 
	);

  //internal function used by align() to update the node path information (for UPLEFT link back nodes ONLY)
  //
  //the point is that we check the upleft node for pathes.
  //if the upleft node is ZERO, we could make a new path and add path to this current one. this is simply
  //if the upleft node is not zero, then we will check all the pathes in the upleft node
  //and we need to compare with this current value, if the current value is the better, we need 
  //create and add the path entry to the current node. we also updated the path information to
  //the vector.
  //also we assume 
  //i is the pattern and j is the subject , current score.
  //lenP and lenS is the length of pattern and subject respcectively.
  void updateNodePath(const unsigned& i, const unsigned& j, const unsigned& lenP, 
		const unsigned& lenS, const double & compval);

	//this is called each time we trace back one path. The rationale is that after the tracing back of 
	//of one path (current_pathID), we have used some nodes in the table. these nodes can have 
	//multiple links, so in this way, some other path become dead and not valid anymore. For 
	//those pathes that were involved we have "remember" them in a vector (vec_pathToBeChecked).
	//now we need to go through the path vector to update in case 
	//some node has been take by the current path. this is the driving function will call
	//the recursive function to update the whole path for those remebered during the tracing back. 
	//input: 
	//	current_pathID, the ID for the currently traced back.
	//	vec_pathToBeChecked, the pathes that has been "remembered" in the previous tracing back.
	//					//
	void updatePathForMultipleLinks(const unsigned& current_pathID, 
				vector<unsigned>& vec_pathToBeChecked);
	
	//this function is used to check whether a path (denoted by a PathID) goes throught
//a node (denoted by index_p and index_s). A path not going through a node because it 
//is go through by other path uniquely with higher score during align() and it could also be caused due to 
//it is a node with multiple pathes, but were taken by previous passing of updating. 
// PathID is the same as the location/index of the path in the original c_path_vec formed after align();
bool isPathThruThisNode(const unsigned& pathID, 
			const unsigned& index_p, const unsigned& index_s) const; 


//function to remove the pathID form the pathelementEntry table pathID_vec for the specific node
//   note, so that the node not in the path anymore, mainly because the nodes has been taken and no more
//available now.

void removePathFromNode(const unsigned& pathID, 
					const unsigned& index_p, const unsigned& index_s);

//
//----Recursive--- function
//input: removThisPath, boolean, used to indicate that we need to remove the this path from this node
//		and also all the nodes following this current node, since the road is cut off 
//		anyway.
//	PathID, this current path being checked.
//	optimal score is the best score so far for the path being checked.  
void updatePathForMultipleLinks_Node(const unsigned& pathID, const unsigned& index_p, const unsigned& index_s
				, const double& optimalScore, const bool& removeThisPath, const LinkBack& link);


  //**********************************
  //the following are the ones used to return and keep track of all non-intersect alignments
  //this is the array holding numOfAlignments requested
  //the one defined by base class only holds the optimal one
  unsigned int c_numOfAlignments;
  AlignmentString* c_alignmentArr; //in this array, we hold the best a few of alignments. the number of best ones are defined by input.
  //Note: c_alignment holds the best alignment.
  double* c_scoreArr;
  //in this class we define the alignmentstring and score string and
  //then delete it upon destruction. so the outside caller need to take care(copy)
  //if they need to use the score or alignment after the alignment scope expires.

  vector<Path*> c_path_vec; //information about path. the order after they created in
	//align() is the path ID. 
  vector<unsigned> c_pathID_vec; //this is added 1/23/2020 and used to remember the pathIDs,
			//to be sorted etc; the deal is that we will keep c_path_vec unchanged after
			//align() and then we will only work on c_pathID_vec sorted and updated
			//since we want to keep path vector in their original order. pathID is
			//only the original index in the c_path_vec. this way we could quickly
			//access it. 
  void traceBackMultiple();//doing trace back to return multiple local alignment
  PathElementEntry* c_PathElementTable;  //hold the information for path element for each nodes.
							//has same layout/structure as the tracetable, so we need this 
							//for the purposes of the doing trace back for the special 
							//cases where we might have multiple linkback/path with gap penalty and
							//mismatch score identical.
  unsigned int c_current_path_tk; //the current path that is checking by traceback()function
  
  
};


#endif

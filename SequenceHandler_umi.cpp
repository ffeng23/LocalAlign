#include <algorithm>    // std::copy
#include <vector>       // std::vector

#include "SequenceHandler.hpp"
#include "Accessory/FastaHandler.hpp"
#include "Accessory/FastqHandler.hpp"
#include "Accessory/string_ext.hpp"
#include "OverlapAlignment.hpp"
#include "GlobalAlignment.hpp"
#include "LocalAlignment.hpp"
#include "SequenceHandler_umi.hpp"


//in this function, we pass in an alignment string and get the alignment strings with gapextension
//then cout and extract the umi string from it.
//we check for each subject (bottom string) for Ns and get the appropriate UMI from the top string.
//it could be '-'.
//static 
//----- OBSLETED--------
string GetUmiFromAlignmentString(const AlignmentString& aso)
{
	string umi;
	string p_w_gap;
	string s_w_gap;
	p_w_gap=aso.GetPattern(true);
	s_w_gap=aso.GetSubject(true);
	
	for(unsigned i=0;i<s_w_gap.size();i++)
	{
		char s=s_w_gap.at(i);
		if(s=='N')
		{
			umi.append(1,p_w_gap.at(i));
		}
	}
	//return(move(umi));
	return umi;
}

//based on the mapped umi_anchor position, trim the sequence from the beginning. 
//input umi_anchor-end is the mapped length based on the umi and anchor.
//static
 string trimSequenceUmiAnchorFromStart(const string& seq, const unsigned &umi_anchor_end)
{
	if(umi_anchor_end+1>seq.size())
	{
		return "";
	}
	return (move(seq.substr( umi_anchor_end+1)));
}

//we basically extract positions of umi. 
//for now we extend the vocabulary to have N and H, both of which are treated as umi.
//others are treated as spacer or anchor
//input: 
//          string umi_pattern character string
//Output:
//              umi_positions, umi_length
//              anchor_positions, anchor_length 
bool parseUmiPattern(const  string &umi_pattern, 
                    unsigned* const umi_positions, unsigned &umi_length,
					 unsigned* const anchor_positions, unsigned &anchor_length)
{
	unsigned plen=umi_pattern.size();
	umi_length=0;
	anchor_length=0;
	//check each individual char in the pattern
	for(unsigned i=0;i<plen;i++)
	{
		char c=umi_pattern.at(i);
		if(!is_nucleotide(c))
		{
			cerr<<"ERROR: contains invalid nucleotide code\'"<<c<<"\' at position "<<i<<endl;
			return false;
		}
		if(c=='N'||c=='n'||c=='H'||c=='h')
		{
			umi_positions[umi_length]=i;
			umi_length++;
		}
		else
		{
			anchor_positions[anchor_length]=i;
			anchor_length++;
		}
	}
	if(umi_length>umi_pattern.size()||anchor_length>umi_pattern.size())
	{
		cerr<<"ERROR: the umi length or anchor lenght is larger then umi_pattern string, invalid1!!"<<endl;
		return false;
	}
	return true;				 
}
//follow umi_tools style to insert the umi barcode to the sequence name. 
const string insertUmi2SName_umitools(const string& sname, const string& umi)
{
	unsigned in=sname.find_first_of(' ');
	if(in==std::string::npos)
	{
		return sname;
	}
	string temp(sname.substr(0,in));
	temp.append("_T_");
	temp.append(umi);
	temp.append(sname.substr(in));
	//sname.insert(in, umi.insert(0,"_T_"));
	//return move(temp);
	return temp;
}
//we prepare the anchor string for aligning.
//we assume anchor will not be zerolength
//we strip off the Ns at the both end of the anchor. and leave the 
//middle Ns in the anchor.  <--- modified now (8/28/21) in order to take care cases: NNGGNNGGG
//          in which there are gaps between spacer.  So now we on the basis of anchor positions to pick the exact space and
//              get rid of Ns in the middle.
//anchor_length is actuall number of non N code. it is not the anchor span.
const string getAnchor(const string &umi_pattern, 
			const unsigned* const anchor_positions, const unsigned & anchor_length)
{
	if(anchor_length==0)
		return "";
	string anchor ;
	//get rid of the right extra Ns
	//anchor.assign(move(anchor.substr(0, anchor_positions[anchor_length-1]+1)));
	//anchor.assign(move(anchor.substr(anchor_positions[0])));
    for(unsigned i =0;i<anchor_length;i++)
    {
           anchor.push_back(umi_pattern.at(anchor_positions[i]))        ;
    }
	//return move(anchor);
	return move(anchor);
}
/*
//we prepare the umi string for insertion.
//we assume umi will not be zerolength
//we strip off the anchor and leave the 
//Ns in the UMI.
//
const string getUMI(const string &umi_pattern, 
			const unsigned* const umi_positions, const unsigned & umi_length)
{
	string umi;
	for(unsigned i=0;i<umi_length;i++)
	{
		umi.append(1,umi_pattern.at(umi_positions[i]));
	}
	return move(umi);
}*/

//modified 8/28/2021, doing prepare aligned sequence by each segment.
//          is the part of 
//To get the alignment string (part of the squence we want to get umi out of ). we add extra nucletide (by a number of offset) to both ends.
//          now (8/28/2021) we have to expand this, we need to take care of the case where there are gaps in between 
//          the spacers.   NNATTNNHNNATTNNACATGGG. In this case, we will have to get all combinations of offsetted alignSequences
//              no just both ends, but also the gaps. for both ends we simply add extra as we originally did.
//          To take care of the gaps we need to do all combinations of offset =1 at each edge of the spacers.
//                  eg.    NNATTNN...
//                  to take care of ATT at its right edge, we will do offset of +/-1
//              so we will have ATT, ATT(N), AT, ...
//                  in the end, we will have an array of sequences that will be aligned against by spacer/anchor sequences.
//              The good news is that since we have consider all combinations, we only need to do matching against spacers and then
//              and pick the best one. We don't have to do alignment.  (update 8/30/2021: NO: we still do alignment,
//                      in this case, we allow in/del on spacers,. the above one is actually taking care of in/del in gaps.
//                              we will do overlap alignment!!! (ignore the end for score.); in addition this will take care of
//                                  off offset shifting for matching.like indel in the beginning of spacer.)
//            To realize above we actually work on the gap segment( on input sequence not umi), ATT only place holder indicating 
//                  spacer locations.
//                    ATT__ATT   => perfect case: ATTATT  + NN (UMI), this case is not distinguish from AT___ATT
//                    ATT___ATT   =>ATTATT + NNN (UMI), not distiguishable from AT____ATT 
//                    ATT_ATT  =>ATT_ATT + N (UMI)   this one is not distiguishable from ATT_XTT  or ATX_ATT
//     again, it is hard to distinguish deletion from mutation, for example a deletion in first ATT lead to 
//                                  ATT__ATT --> AT__ATT,  which could be detected as  ATX_ATT in the case III with a mutation.  
//  we assume there will be in/del in spacer sequences.
//           On each  gap we will allow one offset. This is only one gap. this will combine with other gap (full combinatorial)
//but there are special cases, where we are dealing there are no more extra to add on the end. 
//since adding N on it won't give us a good match.
//we will do overlap alignment, and take care of the case, where we shift toward front/left.
// the special cases are dectected as 
// 				anchor_positions[0]-offset<0
//				anchor_positions[last]+offset>=seq.length() 
void prepAlignSequence(const string& seq, const unsigned& offset, //const unsigned& offset_right, 
		const unsigned* const anchor_positions, const unsigned& anchor_length, 
        //const unsigned* gap_start,  
        vector<string>* alignSeq, /*output*/
		vector<vector<unsigned>>* alignSeq_start /*output*/, vector<vector<unsigned>>* alignSeq_end /*output*/)
{
    //vector<string> alignSeq;
    //cout<<"inside prep function......"<<endl;
	if(anchor_length==0)  //no spacer, we will only take a fixed number of UMIs, NNNN from the sequence.
		return ;
    //first figure out how many gaps, we only need to take care middle gaps ATTNNNATT, 
    //but not edge gaps NNATTNN. (NN is a gap)
    //by going through the anchor_positions 
    unsigned gap_num=0;
    vector<unsigned> gap_start; //remember each gap starting 
    vector<unsigned> gap_end;//remeber each gap ending
     //   cout<<"Start doing loop.........."<<endl;
    for(unsigned i=0;i<anchor_length-1;i++)  //this way, we simply "overlook" the end gap, the NNN at both end of the sequence. 
    {
        if(anchor_positions[i+1]-anchor_positions[i]>1)//this is a gap
        {
            gap_num++;
            gap_start.push_back(anchor_positions[i]+1);
            gap_end.push_back(anchor_positions[i+1]-1);
        }
    }
    
    //take care of left edge gap first, if there is one.
    // well it doesn't matter whether there is a gap or not.
	//take care of extras first 
	unsigned left=anchor_positions[0];
	//unsigned right=anchor_positions[anchor_length-1];
	if(left>=offset)
		left=left-offset;
	else
		left=0;
    unsigned right=anchor_positions[0];
    right=right+offset;
    if(right >=seq.size())
           right=seq.size()-1;
    //cout<<"before calling recursive ....... left:"<<left<<"; offset : "<<offset<<endl;
    //cout<<"gap_num:"<<gap_num<<endl;
    //get substring recursively
    //initialize the space as output of the recursive function 
    unsigned current_gap=0;
    for(unsigned i=left;i<=right;i++)
    {
    getSubSeq_recursive(seq, 
                         gap_start, //
                        gap_end, 
                        gap_num,
                        anchor_positions, 
                        anchor_length, 
                     
                     i ,   //current subseq starting position, will be changing pointing to the new subsequence. 
                    /*input &output*/  
                    current_gap, //current gap index, will be changing too.
                         alignSeq, 
                          alignSeq_start,
                         alignSeq_end);
    }
    /*very end edge gap. need to take care here 
	if(right+offset+1<=seq.size())
		right=right+offset;
	else
		right=seq.size()-1;
	seq_str_start=left;
	seq_str_end=right;
	*/
    //cout<<"Done with prep function"<<endl;
	return;
}

//Check the description in preAlignSequence () about the algorithm!!!
//return the string recursively that has concatenated together all proceeding portion.
//      input: Seq , this is the whole string/sequences, that we are working on to get matched and then umis.
//                  
void getSubSeq_recursive(const string& seq, 
                        const vector<unsigned>& gap_start, //
                        const vector<unsigned>& gap_end, 
                        const unsigned & gap_length,
                        const unsigned* const anchor_positions, 
                        const unsigned& anchor_length, 
                     
                      const unsigned& start ,   //current subseq starting position, will not be changing pointing to the new subsequence. 
                     /*input &output*/ 
                     const unsigned& current_gap, //current gap index, will be changing too.
                         vector<string>* alignSeq, 
                         vector<vector<unsigned>>* alignSeq_start,
                            vector<vector<unsigned>>* alignSeq_end
            )
{
        //3 cases: no shift
        //                  left shift  minus one offset
        //                  right shift    plus one offset
        //all recusively
        
            unsigned seq_end=0;
            
            string temp;
           // cout<<"\tcalling resusieve"<<endl;
            //take care special case first.
            //last last gap, we will return now. but before we return we will need to take care stuff first.
            //seqEnd and spacer string. we will also stop on spacer.
            if(current_gap==gap_length)   //including the gap_length==0 case, 
            {
                //cout<<"in final before return:"<<endl;
                if(gap_length==0) //no gap at all only one spacer chunk.
                {
                    seq_end=start+anchor_positions[anchor_length-1]-anchor_positions[0];
                }
                else{ //there are gaps, but we are doing in last one. 
                    seq_end=start+(anchor_positions[anchor_length-1]-gap_end[current_gap-1])-1; //here note the gap_end and anchor_positions are all expected position, could be used for figuring out the relative positions
                }
                temp=seq.substr(start, seq_end-start+1);
                alignSeq->push_back(move(temp));
                vector<unsigned> tempStart;
                tempStart.push_back(start);
                (*alignSeq_start).push_back(tempStart);
                vector<unsigned> tempEnd;
                tempEnd.push_back(seq_end);
                (*alignSeq_end).push_back(tempEnd);
                
                return ;
            }
            
            //cout<<"\t\t not done yet"<<endl;
            //doing the job to call on three different conditions.
            //note: everything is relative. otherwise it is not what we want.
            //case 1, no shifting, the best scenario
            //first we need to get the preceeding spacer length,
            //to do that we check the current gap start and previous gap end. or 
            // or if this is the first gap, then we just take gap start -1 since there is no previous gap and no shifting gaps
            //cout<<"\t\t current_gap:"<<current_gap<<endl;
            if(current_gap==0)   //note, if we are here, means the gap_length>0, since gap_length==0 is consider in the previous condiiton (up)
                                //also it is not possible for anchor_length==0 in here.
            {
                seq_end=start+gap_start[current_gap]-anchor_positions[0]-1;
            }
            else //could be last one. or some other ones not first though. doing relative to the start point
            {
                seq_end=gap_start[current_gap]-gap_end[current_gap-1]-1;
                seq_end+= start-1;
            }
            //cout<<"\t\t start:"<<start<<"; seq_end:"<<seq_end<<endl;
            //get current sub string
            temp=seq.substr(start, seq_end-start+1);
            //cout<<"\t\tcurrent sub seq string:"<<temp<<endl;
            //get thing ready for next recursive call.
            //we are here meaning there could be more gaps following
            
             //case 1, no shifting
            unsigned next_start =seq_end+(gap_end[current_gap]-gap_start[current_gap]+1)+1;
            unsigned next_gap=current_gap+1;
             //cout<<"\t\t next gap:"<<next_gap<<endl;
            vector<string> tempAlignSeq1;
            vector<vector<unsigned>> tempAlignSeq_start1;
            vector<vector<unsigned> > tempAlignSeq_end1;
            //cout<<"\t\tstart calling recursive for case 1, start:"<<next_start<<endl;
            getSubSeq_recursive(seq, gap_start, gap_end,gap_length,
                                    anchor_positions, anchor_length,
                                    next_start, next_gap, 
                                &tempAlignSeq1, &tempAlignSeq_start1, &tempAlignSeq_end1);
            //need to add current information.
            for(vector<string>::iterator it = tempAlignSeq1.begin(); it != tempAlignSeq1.end(); ++it)
            {
                it->insert(0,temp);
            }
            for(vector<vector<unsigned>>::iterator it = tempAlignSeq_start1.begin(); it != tempAlignSeq_start1.end(); ++it)
            {
                it->insert(it->begin(),start);
            }
            for(vector<vector<unsigned>>::iterator it = tempAlignSeq_end1.begin(); it != tempAlignSeq_end1.end(); ++it)
            {
                it->insert(it->begin(),seq_end);
            }
            //doing the second case, shift one left 
            //<----need to consider special case.!!!!, DONE
            next_start =next_start-1;
            //current_gap++;  <-have done in the previous section
            vector<string> tempAlignSeq2;
            vector<vector<unsigned>> tempAlignSeq_start2;
            vector<vector<unsigned>> tempAlignSeq_end2;
            getSubSeq_recursive(seq, gap_start, gap_end,gap_length,
                                    anchor_positions, anchor_length,next_start, next_gap, 
                                &tempAlignSeq2, &tempAlignSeq_start2, &tempAlignSeq_end2);
            
            //need to add current information.
            for(vector<string>::iterator it = tempAlignSeq2.begin(); it != tempAlignSeq2.end(); ++it)
            {
                it->insert(0,temp);
            }
            for(vector<vector<unsigned>>::iterator it = tempAlignSeq_start2.begin(); it != tempAlignSeq_start2.end(); ++it)
            {
                it->insert(it->begin(),start);
            }
            for(vector<vector<unsigned>>::iterator it = tempAlignSeq_end2.begin(); it != tempAlignSeq_end2.end(); ++it)
            {
                it->insert(it->begin(),seq_end);
            }
            
            //doing the third case, shift one right
            //<----- need to consider special case
            next_start=next_start+1+1;
            
            //current_gap++;  <-have done in the previous section
            vector<string> tempAlignSeq3;
            vector<vector<unsigned>> tempAlignSeq_start3;
            vector<vector<unsigned>> tempAlignSeq_end3;
            getSubSeq_recursive(seq, gap_start, gap_end,gap_length, 
                                 anchor_positions, anchor_length,next_start, next_gap,
                                &tempAlignSeq3, &tempAlignSeq_start3, &tempAlignSeq_end3);
            
            //need to add current information.
            for(vector<string>::iterator it = tempAlignSeq3.begin(); it != tempAlignSeq3.end(); ++it)
            {
                it->insert(0,temp);
            }
            for(vector<vector<unsigned>>::iterator it = tempAlignSeq_start3.begin(); it != tempAlignSeq_start3.end(); ++it)
            {
                it->insert(it->begin(),start);
            }
            for(vector<vector<unsigned>>::iterator it = tempAlignSeq_end3.begin(); it != tempAlignSeq_end3.end(); ++it)
            {
                it->insert(it->begin(),seq_end);
            }
            
            //concatenate 3 vectors and return
            std::copy(tempAlignSeq1.begin(), tempAlignSeq1.end(), back_inserter(*alignSeq));
            std::copy(tempAlignSeq2.begin(), tempAlignSeq2.end(), back_inserter(*alignSeq));
            std::copy(tempAlignSeq3.begin(), tempAlignSeq3.end(), back_inserter(*alignSeq));
            
            std::copy(tempAlignSeq_start1.begin(), tempAlignSeq_start1.end(), back_inserter(*alignSeq_start));
            std::copy(tempAlignSeq_start2.begin(), tempAlignSeq_start2.end(), back_inserter(*alignSeq_start));
            std::copy(tempAlignSeq_start3.begin(), tempAlignSeq_start3.end(), back_inserter(*alignSeq_start));
            
            std::copy(tempAlignSeq_end1.begin(), tempAlignSeq_end1.end(), back_inserter(*alignSeq_end));
            std::copy(tempAlignSeq_end2.begin(), tempAlignSeq_end2.end(), back_inserter(*alignSeq_end));
            std::copy(tempAlignSeq_end3.begin(), tempAlignSeq_end3.end(), back_inserter(*alignSeq_end));
            return ;
}


//extract umi for each individual sequence
bool extractUmi(SequenceString & ss, const string& umi_pattern, 
				const unsigned* const anchor_positions, const unsigned& anchor_length,
				const unsigned* const umi_positions, const unsigned& umi_length,
				const bool& trim, const bool& extract, 
				ScoreMatrix*  sm, const double &gapopen, const double &gapextension,
				const double & scale,
				const unsigned & Mismatches, const unsigned & OffsetForward,
				string& ex_umi
				)
{
	//assume we either need to extract or trim or both.
	
	//prepare the umi pattern to make the anchor ready
	//string ex_umi;
	if(anchor_length==0)//no anchor,all umis. this is easy, we simply extract the umi 
	{
		ex_umi.append(ss.GetSequence().substr(0,umi_length));
		if(trim)
			ss.SetSequence(ss.GetSequence().substr(umi_length));
		if(extract)
		{
			ss.SetName(insertUmi2SName_umitools(ss.GetName(), ex_umi));
		}
		return true;
	}
	string anchor;
	if(umi_length==0) //no umi, only have anchor, we can only do trimming, but no extract
	{
		ex_umi.assign("");
		anchor.assign(umi_pattern);
	}
	else
	{
		//umi first
		anchor=getAnchor(umi_pattern, anchor_positions, anchor_length);
	}
	
	//prepare the ss so that we can align to it.
	//vector<unsigned> seq_str_start;
    //vector<unsinged>seq_str_end;
    vector<string> alignSeq;
    vector<vector<unsigned> > alignSeq_start; //used to remember after offsetted/shifted, where the sequence string to start and end 
										//these variable is used so as to know where to search the umi strings.
  vector<vector<unsigned> > alignSeq_end;
  prepAlignSequence(ss.GetSequence(), OffsetForward,  
		anchor_positions, anchor_length,   
        &alignSeq, /*output*/ 
        &alignSeq_start /*output*/, 
        &alignSeq_end /*output*/);
    
    //now we got the aligment string. need to pick the best one by alignment,
    //again, we have considered all combinations of alignment string, now we 
    //only need to do overlap alignment to get the best one.
	//string as(move(prepAlignSequence(ss.GetSequence(), OffsetForward, anchor_positions, 
	//			anchor_length, seq_str_start, seq_str_end)));
    unsigned best_index=alignSeq.size()+1; //best_index starting as one pointing outside the array, meaning no good.
    AlignmentString aso;
    AlignmentString best_aso;
    double score, max_score=5*anchor.size();
    double best_score=-1000;
    cout<<"Total alignSeq size:"<<alignSeq.size()<<endl;
    cout<<"inside extract Umi, gapopen:"<<gapopen<<endl;
    for(unsigned j=0;j<alignSeq.size();j++)
    {
        cout<<"Doing alignment at round :"<<j<<endl;
        SequenceString as_ss("s",alignSeq.at(j));
        SequenceString anchor_ss("anchor", anchor);
        
        /*
            //shifted/offsetted boundary, are similar to above seq_str_start and seq_str_end
        //the difference is that the boundary are the real boundary, could be negative
        //the seq_str_start and _end are the practical one, can not go negative.
        int anchor_offsetted_boundary_start; unsigned anchor_offsetted_boundary_end;
        
        //let get the boundary first
        //note the anchor_offsetted_boundary has nothing to do with alignment
        //we originally have the anchor to end align with the sequence 
        //this is the "ideal" positoin. So the shift is relative this "idea' positoin.
        //After we do the below alignment, we check the best aligned (or two) pair 
        //to check for the shift and offset.
        anchor_offsetted_boundary_start=anchor_positions[0]-OffsetForward;
        anchor_offsetted_boundary_end=anchor_positions[anchor_length-1]+OffsetForward;
        */

        //now do alignment between as and anchor string so to find the location
        OverlapAlignment ol(&as_ss, &anchor_ss, sm, gapopen, gapextension
                            ,scale, 0/*affine gap model*/);
        //GlobalAlignment ol(&as_ss, &anchor_ss, sm, gapopen, gapextension
        //				,scale);
        //unsigned num_alignment=2;
        //LocalAlignment ol(&as_ss, &anchor_ss, sm, gapopen, gapextension
        //				,scale, num_alignment, 0); //calling specifically the affine gapmodel and return 2 best alignment
        aso=ol.GetAlignment(); //pattern is on top; and subject is on bottom.
        //for(unsigned j=0;j<num_alignment;j++)
                    //first get the error/mismatch numbers, we can only get errors from the alignment string (case 2). 
            //no errors due to the stripped off 'N' at the end (case 1, and 3)
            score=aso.GetScore();
            //max_score=5*aso.GetSubject().size();//anchor.size();
            //double errors =(max_score-score)/(4+5);
            cout<<"aso alignment:"<<aso.toString()<<endl;
            cout<<"score:"<<score<<"; max_score:"<<max_score<<endl;
            if(max_score ==score)
            {
                //we are doing so, we need to stop and go with this best one.
               best_score=score;
                best_index=j;
                best_aso=aso;
                break;
            }
            if(score >best_score)
            {
                cout<<"set the best one!!!"<<endl;
                best_score=score;
                best_index=j;
                best_aso=aso;
            }
    }
    if(best_index >alignSeq.size())
    {
        //no good!!!
        ex_umi.assign("");
        return false;
    }
    
    //first get the error/mismatch numbers, we can only get errors from the alignment string (case 2). 
		//no errors due to the stripped off 'N' at the end (case 1, and 3)
		score=best_aso.GetScore();
		//max_score=5*aso.GetSubject().size();//anchor.size();
		double errors =(max_score-score)/(4+5); //5 is match, then a mismatch is -4, so 5+4 is the difference for one mismatch.
		
		cout<<"2score:"<<score<<"; max_score:"<<max_score<<endl;
		//cout<<"error:"<<errors<<endl;
		if(((double)Mismatches)<errors)
		{
			//bad!!!
			ex_umi.assign("");
			//cout<<"\t\ttoo many errors; round :"<<j<<endl;
			//continue;
			return false; 
		}
		//AlignmentString aso=aso_arr[j];
		unsigned p_start=best_aso.GetPatternIndexStart();
		//unsigned p_end=aso.GetPatternIndexEnd();
		unsigned s_start=best_aso.GetSubjectIndexStart();
		//unsigned s_end=aso.GetSubjectIndexEnd();
		cout<<"doing.offset now........."<<endl;
        //consider offset now.  cases are not good if subject not starting at zero, or patter started at 2*offset 
        if(s_start!=0||p_start >2*OffsetForward)
        {
               ex_umi.assign("");
            return false;
        }
        
        //We are here means we are going on get an alignment with pattern.
		//start working through the pattern string to get umis out.
        //first need to be ready for doing extract and/or trim, or both. (but not neither.)
        //we first get ready to know the length of the aligned anchor string on pattern (top),
        // in case that we need to do extract
        cout<<"Checking trim extract............alignSeq_end.szie():"<<alignSeq_end.size()<<endl;
        cout<<"best_index:"<<best_index<<endl;
        unsigned p_end_real=alignSeq_end.at(best_index).at(alignSeq_end.at(best_index).size()-1);  //note: at least one anchor is here, if no anchor we will not be here
        
        if(trim )
        {
            if(!extract )
            {//trim only, we don't have to do extract.
                    ex_umi.assign(""); //still set it as "", even though it is not empty
                    cout<<"Here to extract....."<<endl;
                    ss.SetSequence(trimSequenceUmiAnchorFromStart(ss.GetSequence(), p_end_real));
                    return true;
            }
        }
        cout<<"Get umi...........best_index:"<<best_index<<endl;
        //we have we either need to extract alone or trim and extract 
        //Get UMI by alignment and alignSeq_start and anignSeq_end
        //well it is easier now, we don't have to know about alignment, only checking the
        //alignSeq_start and end, plus the end gaps if they are available.
        //since we have already found the best alignment.
        //the hardest part the end gap or edge umis.
        //first do starting umi if they are available.
        //is there starting umi?? check umi_positions and anchor_positions
        string umi_string;
        string seqString=ss.GetSequence();
        unsigned temp_start;
        unsigned temp_end;
        vector<unsigned> best_alignSeq_start=alignSeq_start.at(best_index);
        vector<unsigned> best_alignSeq_end=alignSeq_end.at(best_index);
        cout<<"best align start size:"<<best_alignSeq_start.size()<<"; best_alignSeq_end size:"<<best_alignSeq_end.size()<<endl;
        cout<<"SeqString:"<<seqString<<endl;
        
        unsigned pattern_align_real_start=0; //this is used to determine, where the pattern aligned part starts. 
                                                                                //this is necessary, since we might have shift. Shift make this is not 
                                                                                //agree with best_alignSeq_start[0] (note this include biggest shift), 
                                                                                //another thing make it more complicated is that it is possible that the
                                                                                //alignement could start from middle (not from the first anchor nt.!!! we need to 
                                                                                //figure out this "shift" too)
        if(umi_length>0)
        {
            if(anchor_length>0)
            {
              //compare them. one of them should be zero 
                if(umi_positions[0]<anchor_positions[0])
                {
                    //figure the substring positions
                    cout<<"anchor_positions[0]:"<<anchor_positions[0]<<"; umi_positions[0]"<<
                                    "best_alignSeq_start[0]:"<<best_alignSeq_start[0]<<
                                    "best_alignSeq_end[0]:"<<best_alignSeq_end[0]<<endl;
                    cout<<"best_aso.GetPatternIndexStart():"<<best_aso.GetPatternIndexStart()
                            <<"; best_aso.GetSubjectIndexStart():"<<best_aso.GetSubjectIndexStart()<<endl;
                    //first determine the real pattern aligned starting point. (Note: 9/5/2021, this is too compliated. might not necessary, since we have modified the preAlignSequence part. but this is still right.)
                    int fromStart=best_aso.GetPatternIndexStart()-best_aso.GetSubjectIndexStart(); //note the indexes are relative
                    pattern_align_real_start=fromStart+best_alignSeq_start[0]; //note best_alignSeq_start are absolute, so now e have absolute.
                                    //also very importantly, this is starting of the subject, not the aligned part. and this is what we want!!!
                    if(((int)pattern_align_real_start)<0)
                        pattern_align_real_start=0;
                    temp_start=anchor_positions[0]-umi_positions[0]; //relative length for umi. how long the starting edge umi should be
                    //temp_start=0;
                    //in case there are shift, we need to take care of this
                    if(pattern_align_real_start>=temp_start)
                    {
                        temp_start=pattern_align_real_start- temp_start;
                    }
                    else
                    {
                        temp_start=0;
                    }
                    temp_end=0;
                    if(pattern_align_real_start>1)
                        temp_end=pattern_align_real_start-1;
                    cout<<"temp_start:"<<temp_start<<"; temp_end:"<<temp_end<<endl;
                    umi_string.append(seqString.substr(temp_start, temp_end-temp_start+1));
                    cout<<"starting edge umi:"<<umi_string<<endl;
                }
                else  //in this case, there is no edge umis, it starts with anchor.
                {
                    //do nothing.
                }
            }
            else //can not happen 
            {
                cerr<<"something wrong, can not be both zero."<<endl;
                exit(-1);
            }
        }
        else // umi_length zero, can not happen.   
        {
            cerr<<"something wrong, can not be umi length zero."<<endl;
            if(anchor_length<=0)
            {
                cerr<<"something wrong, can not be both zero."<<endl;
                exit(-1);
            }
            exit(-1);
        }
        cout<<"doing middle umi.............."<<endl;
        //now doing the middle umis if there are any.
        for(unsigned k=0;k<best_alignSeq_start.size()-1;k++)
        {
            cout<<"Looping throught middle pieces...."<<endl;
            temp_start=best_alignSeq_end.at(k)+1;
            temp_end=best_alignSeq_start.at(k+1)-1;
            umi_string.append(seqString.substr(temp_start, temp_end-temp_start+1));
        }
        
        cout<<"umi string.:"<<umi_string<<endl;
        cout<<"doing end umi.............."<<endl;
        //now doing the last piece if there is one.
        if(umi_positions[umi_length-1]>anchor_positions[anchor_length-1])
		{
            temp_start=best_alignSeq_end.at(best_alignSeq_end.size()-1)+1;
            
            //take care of the special case where temp_end is too long due to shift.
            temp_end=umi_positions[umi_length-1]-anchor_positions[anchor_length-1]; // realtive!!!!
            if(temp_start+temp_end >seqString.size() ) //again too long due to shifting.
            {
                temp_end=seqString.size()-temp_start;
            }
            umi_string.append(seqString.substr(temp_start, temp_end));
            
            //also need to update p_read_end, in case we need to do trimming
            p_end_real+=temp_end;
        }

		ex_umi.assign(move(umi_string));
		//ex_umi.append(move(umi_inside));
		//ex_umi.append(move(umi_after));
		
		//now let's take care of insertion and and trim
		if(extract)
		{
			ss.SetName(insertUmi2SName_umitools(ss.GetName(), ex_umi));
		}
		if(trim)
		{
            
			//use the previously caculated umi_anchor_end to cat off the sequence
			ss.SetSequence(trimSequenceUmiAnchorFromStart(ss.GetSequence(), p_end_real));
		}
		//if we are here, we are good. 
		return true;
	}	


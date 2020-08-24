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

//we basically
bool parseUmiPattern(const  string &umi_pattern, unsigned* const umi_positions, unsigned &umi_length,
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
		if(c=='N'||c=='n')
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
//middle Ns in the anchor.
//anchor_length is actuall number of non N code. it is not the anchor span.
const string getAnchor(const string &umi_pattern, 
			const unsigned* const anchor_positions, const unsigned & anchor_length)
{
	if(anchor_length==0)
		return "";
	string anchor (umi_pattern);
	//get rid of the right extra Ns
	anchor.assign(move(anchor.substr(0, anchor_positions[anchor_length-1]+1)));
	anchor.assign(move(anchor.substr(anchor_positions[0])));
	//return move(anchor);
	return anchor;
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
//we for the alignment string (top, first). we add extra nucletide (by a number of offset) to both ends.
//but there are special case, where we are dealing there are no more extra to add on the end. 
//since adding N on it won't give us a good match.
//we will do overlap alignment, and take care of the case, where we shift toward front/left.
// the special cases are dectected as 
// 				anchor_positions[0]-offset<0
//				anchor_positions[last]+offset>=seq.length() 
const string prepAlignSequence(const string& seq, const unsigned& offset, 
		const unsigned* const anchor_positions, const unsigned& anchor_length,
		unsigned& seq_str_start /*output*/, unsigned& seq_str_end /*output*/)
{
	if(anchor_length==0)
		return "";
	//take care of extras first 
	unsigned left=anchor_positions[0];
	unsigned right=anchor_positions[anchor_length-1];
	if(left>=offset)
		left=left-offset;
	else
		left=0;
	if(right+offset+1<=seq.size())
		right=right+offset;
	else
		right=seq.size()-1;
	seq_str_start=left;
	seq_str_end=right;
	//cout<<"left: "<<seq_str_start<<"--right:"<<seq_str_end<<endl;
	return(seq.substr(left, right-left+1));
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
	unsigned seq_str_start, seq_str_end; //used to remember after offsetted/shifted, where the sequence string to start and end 
										//these variable is used so as to know where to search the umi strings.
	string as(move(prepAlignSequence(ss.GetSequence(), OffsetForward, anchor_positions, 
				anchor_length, seq_str_start, seq_str_end)));
	SequenceString as_ss("s",as);
	SequenceString anchor_ss("anchor", anchor);
	
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
	

	//now do alignment between as and anchor string so to find the location
	//OverlapAlignment ol(&as_ss, &anchor_ss, sm, gapopen, gapextension
	//					,scale);
	//GlobalAlignment ol(&as_ss, &anchor_ss, sm, gapopen, gapextension
	//				,scale);
	unsigned num_alignment=2;
	LocalAlignment ol(&as_ss, &anchor_ss, sm, gapopen, gapextension
					,scale, num_alignment, 0); //calling specifically the affine gapmodel and return 2 best alignment
	AlignmentString* aso_arr=ol.GetAlignmentArr(); //pattern is on top; and subject is on bottom.
	for(unsigned j=0;j<num_alignment;j++)
	{
		AlignmentString aso=aso_arr[j];
		unsigned p_start=aso.GetPatternIndexStart();
		unsigned p_end=aso.GetPatternIndexEnd();
		unsigned s_start=aso.GetSubjectIndexStart();
		unsigned s_end=aso.GetSubjectIndexEnd();
		cout<<"aso alignment:"<<aso.toString()<<endl;
		//start working through the pattern string 
		//the real pattern start at, the orignal output of alignment is relative to the aligned/searchable sequence string 
		p_start=p_start+seq_str_start;
		p_end=p_end+seq_str_start;
	/*  //This fowllowing is still true, but we are not going to use this.
		//we can only have 3 cases, 1) subject in fall in middle of pattern, this is the easiest case
		//2)subject is left cut;3)subject is right cut.4)sujbect is both cut. this is possible event_callback
		//subject is shorter than pattern; but we can fill gaps in subject.
		if(s_start==0)
		{
			if(s_end==anchor.length()-1) //easies case, full length anchor
			{   
				//now get the sequence string ready
				//Get actually position for pattern aligned
				p_start=p_start+seq_str_start;
				
			}
			else  //in this case, s_end can only smaller than anchor length 
			{ 
			}
		}
		else  //s_start>0, so means on the left, some nt of subject was not aligned/cut off 
		{
			if(s_end<anchor.length()
			{
			}
		}*/
		
	//this is another "new" way, we "calculate" based on offset and alignment 
	// where to start getting the umi, then we "count" till the end and collect
	// the umi and other information, offset/shift/mismatches. and make descision
	//about whether to accept the umi on the fly.
		//Note: we declar the start to be signed int, because it can go negative because
		//shifting and aligning.
		//int umi_anchor_start;   <-now we don't use it. since it doesn't matter. we will not use it.
		unsigned umi_anchor_end;  //these are for the full umi_anchor string, we will calculate how this mapped to the sequence string
		int anchor_start; unsigned anchor_end;
		//unsigned aligned_subject_start; unsigned aligned_subject_end;
			//the above are for the subject(umi_anchor string). we will calculat them
			//based on the alignment string and map them onto the squence string.
		
		//remember on the sequence string, we have 
		//p_start and p_end are the aligned pattern
		//seq_str_start and seq_str_end are the searchable sequence string (for anchor) based on
		//offset and original anchor string.
		
		
		//now calculate the quantities
		//aligned_subject_start=p_start;//p_start has been updated, not relative to searchable seq string;
		//aligned_subject_end=p_end;//p_end has been updated, not relative to searchable seq string
		
		anchor_start=p_start-s_start;
		anchor_end=anchor_positions[anchor_length-1]-anchor_positions[0]-s_end +p_end;  //based on the alignment string.
		
		//umi_anchor_start=p_start-s_start-anchor_positions[0];
		unsigned temp=anchor_positions[anchor_length-1];
		if(umi_length>0&&umi_positions[umi_length-1]>anchor_positions[anchor_length-1])
		{
			temp=umi_positions[umi_length-1];
		}
		umi_anchor_end=temp - anchor_positions[0]-s_end+p_end;
		cout<<"anchor_start:"<<anchor_start<<"; anchor end:"<<anchor_end<<endl;
		cout<<"shifted anchor_start:"<<anchor_offsetted_boundary_start
			<<"; shfited anchor end:"<<anchor_offsetted_boundary_end<<endl;
		
		//now check for offset. hopefully not shifted two much
		if(anchor_start<anchor_offsetted_boundary_start)
		{
			cout<<"shifted two mcuh, round:"<<j<<endl;
			continue;
			//return false; //shifted too left
		}
		if(anchor_end>anchor_offsetted_boundary_end)
		{
			cout<<"shifted two mcuh 2, round:"<<j<<endl;
			continue;
			//return false; //shifted too right
		}
		
		//first get the error/mismatch numbers, we can only get errors from the alignment string (case 2). 
		//no errors due to the stripped off 'N' at the end (case 1, and 3)
		double score=aso.GetScore();
		double max_score=5*aso.GetSubject().size();//anchor.size();
		double errors =(max_score-score)/(4+5);
		
		cout<<"score:"<<score<<"; max_score:"<<max_score<<endl;
		//cout<<"error:"<<errors<<endl;
		if(((double)Mismatches)<errors)
		{
			//bad!!!
			ex_umi.assign("");
			cout<<"\t\ttoo many errors; round :"<<j<<endl;
			continue;
			//return false; 
		}
		//now let us do umi extracting
		//if we are here, meaning we did not shifted too much, go on to get umis and error/mismatches
		//we will do 3 cases. 1) before anchor string;2)within anchor string;3)after anchor string
		//umi position is determined by the starting umi_anchor_start position and umi position.
		//first go through the list to pick up the 3 cases for umi
		//unsigned* umi_before_pos=new unsigned [umi_length];
		unsigned umi_before_length=0;
		string umi_before;
		
		//unsigned* umi_after_pos=new unsigned [umi_length];
		unsigned umi_after_length=0;
		string umi_after;
		
		//unsigned* umi_inside_pos=new unsigned [umi_length];
		unsigned umi_inside_length=0;
		string umi_inside;
		
		int curr_umi_pos=0;
		for(unsigned i=0;i<umi_length;i++)
		{
			if(umi_positions[i]<anchor_positions[0])
			{
				//before
				curr_umi_pos=umi_positions[i];
				//turn into mapped position (based on the alignment strings)
				curr_umi_pos =curr_umi_pos-(anchor_positions[0]+s_start)+p_start;
				if(curr_umi_pos<0)
					umi_before.append("N");
				else
					umi_before.append(1,ss.GetSequence().at(curr_umi_pos));
				umi_before_length++;	 
			}
			else
			{
				if(umi_positions[i]>anchor_positions[anchor_length-1])
				{
					//after
					curr_umi_pos=umi_positions[i];
					//turn into mapped position (based on the alignment strings)
					curr_umi_pos =curr_umi_pos-(anchor_positions[0]+s_end)+p_end;
					if(curr_umi_pos>=((int)ss.GetSequence().length()))
						umi_after.append("N");
					else
						umi_after.append(1, ss.GetSequence().at(curr_umi_pos));
					umi_after_length++;	 
				}
				else
				{
					//inside, we will do this outside of this.
					//umi_inside_pos[umi_after_length]=umi_positions[i];
					umi_inside_length++;
					
					continue;
				}
			}
		}
		
		//now we are hear, we have done before and after if there are. 
		//now the inside alignment string
		if(umi_inside_length>0)
		{
			umi_inside.assign(GetUmiFromAlignmentString(aso));
		}
		ex_umi.assign(move(umi_before));
		ex_umi.append(move(umi_inside));
		ex_umi.append(move(umi_after));
		
		//now let's take care of insertion and and trim
		if(extract)
		{
			ss.SetName(insertUmi2SName_umitools(ss.GetName(), ex_umi));
		}
		if(trim)
		{
			//use the previously caculated umi_anchor_end to cat off the sequence
			ss.SetSequence(trimSequenceUmiAnchorFromStart(ss.GetSequence(), umi_anchor_end));
		}
		//if we are here, we are good. 
		return true;
	}	
	//we are now outside of the loop, still did not return true. that means we are no good.
	return false;
}

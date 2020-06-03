#ifndef SEQUENCEHANDLER_UMI_HPP
#define SEQUENCEHANDLER_UMI_HPP

using namespace std;
//used to parse the umi 
bool parseUmiPattern(const string& umi_pattern, unsigned* const umi_positions, unsigned &umi_length,
					 unsigned* const anchor_positions, unsigned &anchor_length);

bool extractUmi(SequenceString & ss, const string& umi_pattern, 
				const unsigned* const anchor_positions, const unsigned& anchor_length,
				const unsigned* const umi_positions, const unsigned& umi_length,
				const bool& trim, const bool& extract, 
				ScoreMatrix* sm, const double &gapopen, const double &gapextension,
				const double & scale,
				const unsigned & Mismatches, const unsigned & OffsetForward,
				string& ex_umi
				);

const string insertUmi2SName_umitools(const string& sname, const string& umi);

//we prepare the anchor string for aligning.
//we assume anchor will not be zerolength
//we strip off the Ns at the both end of the anchor. and leave the 
//middle Ns in the anchor.
//
const string getAnchor(const string &umi_pattern, 
			const unsigned* const anchor_positions, const unsigned & anchor_length);
/*			
//we prepare the umi string for insertion.
//we assume umi will not be zerolength
//we strip off the anchor and leave the 
//Ns in the UMI.
//
const string getUMI(const string &umi_pattern, 
			const unsigned* const umi_positions, const unsigned & umi_length);*/
//we for the alignment string (top, first). we add extra nucletide (by a number of offset) to both ends.
//but there are special case, where we are dealing there are no more extra to add on the end. 
//since adding N on it won't give us a good match.
//we will do overlap alignment, and take care of the case, where we shift toward front/left.
// the special cases are dectected as 
// 				anchor_positions[0]-offset<0
//				anchor_positions[last]+offset>=seq.length() 
const string prepAlignSequence(const string& seq, const unsigned& offset, 
		const unsigned* const anchor_positions, const unsigned& anchor_length,
		unsigned& seq_str_start /*output*/, unsigned& seq_str_end /*output*/
		);
string trimSequenceUmiAnchorFromStart(const string& seq, const unsigned &umi_anchor_end);
string GetUmiFromAlignmentString(const AlignmentString& aso);
#endif
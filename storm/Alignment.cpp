/*
 * Alignment.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "Alignment.h"



Alignment::Alignment()
{
	// TODO Auto-generated constructor stub
	subjectRead = NULL;
	queryRead = NULL;
	queryOrientation = true;
	subjectOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}
Alignment::Alignment(SubjectRead* sRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = sRead;
	queryRead = NULL;
	queryOrientation = true;
	subjectOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}
Alignment::Alignment(QueryRead* qRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = NULL;
	queryRead = qRead;
	queryOrientation = true;
	subjectOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}

Alignment::Alignment(SubjectRead* sRead, QueryRead* qRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = sRead;
	queryRead = qRead;
	queryOrientation = true;
	subjectOrientation = true;
	subjectStart = 0;
	queryEnd = 0;
	subjectEnd = 0;
	editInfor = NULL;
	containedTag =0;

}

Alignment::~Alignment() {
	// TODO Auto-generated destructor stub
	if(editInfor!=NULL)
	{
		(*editInfor).clear();
		delete editInfor;
	}
	queryRead = NULL;
	subjectRead = NULL;
//	if(this->subjectRead!=NULL)
//		delete this->subjectRead;
}

//******OMEGA original definition for the alignment orientation******
//THis is only for the key matching orientation
	// orient 0
	//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
    // omega coordiate:                              0

	// orient 2
	//	 >--------MMMMMMMMMMMMMMM*************> 			read1
	//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
    // omega coordiate:                              0

	// orient 1 (should be 3 in the edge orientation)
	//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
	//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
    //  0       omega coordiate

	// orient 3 (should be 1 in the edge orientation)
	//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
	//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
    //  0       omega coordiate
//EDGE ORIENTATION
// Orientation of overlap edges
// 0 = u<-----------<v
// 1 = u<----------->v
// 2 = u>-----------<v
// 3 = u>----------->v

//*******************************************************************


	// coordinates of the overlap alignment, which is defined by the query reads
	//  The following case will have a positive subject start position
	//  query:        XXXXXXMMMMMMMM
	//	subject:            MMMMMMMMXXXXXXX
	//  The following case will have a negative subject start position
	//  query:        MMMMMMXXXXXXXX
	//	subject: XXXXXMMMMMM



//the query strand can be forward or reversed, since both cases are hashed in the hash table.
int Alignment::orientationTranslate()
{


	if(this->subjectStart>=0&&this->subjectEnd<=this->queryEnd)return -1;
	if(this->subjectStart<=0&&this->subjectEnd>=this->queryEnd)return -1;
	if(this->subjectStart>0)
	{
		if(this->queryOrientation==true)
			{
			if(this->subjectOrientation==true)
				return 3; //read1 = subject, read2 = query
			else
				return 2;
			}
		else
			{
			if(this->subjectOrientation==true)
				return 1; //read1 = subject, read2 = query
			else
				return 0;
			}
	}
	else if(this->subjectStart<0)
	{
		if(this->queryOrientation==true)
			{
			if(this->subjectOrientation==true)
				return 0; //read1 = subject, read2 = query
			else
				return 1;
			}
		else
			{
			if(this->subjectOrientation==true)
				return 2; //read1=subject, read2 = query
			else
				return 3;
			}
	}
	return -1;
}

string Alignment::editInfoToOmegaString()
{
	//example:
	//#74AT:73TA:72GC:71-G:66G-
	//#NA (no errors)
	string output;
	if(this->editInfor==NULL || this->editInfor->size()==0)
//		output = "#NA";
		output = "NA";
	else
	{
		std::stringstream sstm;
//		sstm<<"#";
		for(std::map<int, char>::iterator it=this->editInfor->begin(); it!=this->editInfor->end(); ++it)
		{
			// coordinates are defined by the query reads
			// all insertion and deletion is on the subject read
			// upper case char means substitution
			// lower case char means insertion
			// 'D' means deletion

			char base = it->second; int position = it->first;
			char queryBase;
			if(this->queryOrientation==true)
				queryBase = this->queryRead->getSequence().at(position);
			else
				queryBase = this->queryRead->reverseComplement().at(position);


			if(it!=this->editInfor->begin())
				sstm<<":";
			if(base=='D')
			{
				//remap the position when the alignment needs to be flipped when the original alignment is stretching to the left
				int omegaposition = position;
				if(this->subjectStart<0)
				{
					omegaposition = -position + this->queryEnd;
					queryBase = this->reverseComplement(queryBase);
				}
				sstm<<omegaposition<<queryBase<<"-";
			}
			else if(base=='A'||base=='T'||base=='C'||base=='G')
			{
				int omegaposition = position;
				if(this->subjectStart<0)
				{
					omegaposition = -position + this->queryEnd;
					queryBase = this->reverseComplement(queryBase);
					base = this->reverseComplement(base);
				}
				sstm<<omegaposition<<queryBase<<base;
			}
			else if(base=='a'||base=='t'||base=='c'||base=='g')
			{
				int omegaposition = position;
				base = toupper(base);
				if(this->subjectStart<0)
				{
					omegaposition = -position + this->queryEnd;
					base = this->reverseComplement(base);
				}
				sstm<<omegaposition<<"-"<<base;
			}

		}
		output = sstm.str();
	}
	return output;
}

void Alignment::reverseAlignment()
{
	this->queryOrientation = ! (this->queryOrientation);
	this->subjectOrientation = ! (this->subjectOrientation);

	int substart = this->queryEnd-this->subjectEnd;
	int subend = this->queryEnd-this->subjectStart;

	this->subjectEnd = subend;
	this->subjectStart = substart;

	//there is no need to change queryStart = 0, and queryEnd = querylength -1;
	if(this->editInfor!=NULL && this->editInfor->size()!=0)
	{
		map<int, char>* newEditInfo = new map<int, char>();
		for(std::map<int, char>::iterator it=this->editInfor->begin(); it!=this->editInfor->end(); ++it)
		{
			// coordinates are defined by the query reads
			// all insertion and deletion is on the subject read
			// upper case char means substitution
			// lower case char means insertion
			// 'D' means deletion

			char base = it->second; int position = it->first;
			char newbase=base;
			int newposition = this->queryEnd-position;

			if(base=='D')
			{
				//do nothing
			}
			else if(base=='A'||base=='T'||base=='C'||base=='G')
			{

				newbase = this->reverseComplement(base);

			}
			else if(base=='a'||base=='t'||base=='c'||base=='g')
			{

				base = toupper(base);
				base = this->reverseComplement(base);
				newbase = tolower(base);
			}
			newEditInfo->insert(std::pair<int, char>(newposition, newbase));

		}
		this->editInfor->clear();
		delete this->editInfor;
		this->editInfor = newEditInfo;

	}
}


string Alignment::getSubjectSubSequence(int alignstart, int alignend)
{
	string outstring = "";
	if(alignstart>=this->subjectStart&&alignend<=this->subjectEnd)
	{
		if(this->subjectOrientation == true)
		{
			int start = alignstart-this->subjectStart;
			int end = alignend - this->subjectStart;
			outstring = this->subjectRead->getSequence().substr(start, end-start+1);
		}
		else
		{
			int start = this->subjectEnd - alignend;
			int end = this->subjectEnd - alignstart;
			string temp = this->subjectRead->getSequence().substr(start, end-start+1);
			outstring = QueryRead::reverseComplement(temp);
		}
	}
	return outstring;
}
int Alignment::reverseCoordinate(int oldcoordinate)
{
	return this->queryEnd-oldcoordinate;
}

int Alignment::distanceBetweenEditInfo(Alignment* align)
{
	if(this->editInfor==NULL || this->editInfor->size()==0)
	{
		if(align->editInfor!=NULL)
			return align->editInfor->size();
		else return 0;
	}
	else
	{
		if(align->editInfor==NULL || align->editInfor->size()==0)
			return this->editInfor->size();
		else
		{
			int unionnumber = calculateWenn(this->editInfor,align->editInfor);
			return unionnumber;
		}
	}
}

int Alignment::calculateWenn(map<int, char>* editInfo1, map<int, char>* editInfo2)
{
	int n1, n2, overcount, n;
	n1 = editInfo1->size();
	n2 = editInfo2->size();
	overcount = 0;
	for(std::map<int, char>::iterator it=editInfo1->begin(); it!=editInfo1->end(); ++it)
	{
				// coordinates are defined by the query reads
				// all insertion and deletion is on the subject read
				// upper case char means substitution
				// lower case char means insertion
				// 'D' means deletion

		char base = it->second; int position = it->first;
		map<int,char>::const_iterator it2 = editInfo2->find(position);
		if(it2!=editInfo2->end())
		{
			if(it2->second!=base)
				overcount = overcount+1;
			else
				overcount = overcount+2;
		}
	}
	n = n1+n2-overcount;
	return n;
}

char Alignment::reverseComplement(char base)
{
	char reverseBase;
	if( base & 0X02 ) // C <==> G
		reverseBase = base ^ 0X04;
	else // A <==> T
		reverseBase = base ^ 0X15;
	return reverseBase;
}

//either subject is contained in query or query is contained in subject
bool Alignment::isContainedAlignment()
{
	if(this->subjectStart<=0&&this->subjectEnd>=this->queryEnd)
	{
		containedTag =2;
		return true;
	}
	else if(this->subjectStart>=0&&this->subjectEnd<=this->queryEnd)
	{
		containedTag =1;
		return true;
	}
	else return false;
}

//Either Left Side Alignment or Right Side Alignment, if the query is always on the positive direction, assuming it shouldn't be contained alignment
bool Alignment::isLeftOrRightAlignment(bool& isRight, int& startPoint)
{
	if(this->queryOrientation==true)
	{
		if(this->subjectStart>=0)
		{
			isRight = true;
			startPoint = this->subjectStart;
		}
		else
		{
			isRight = false;
			startPoint = this->queryEnd - this->subjectEnd;
		}

	}
	else
	{
		if(this->subjectStart>=0)
		{
			isRight = false;
			startPoint = this->subjectStart;
		}
		else
		{
			isRight = true;
			startPoint = this->queryEnd - this->subjectEnd;
		}
	}
	return true;
}

int Alignment::getEditDistance()
{
	int size;
	if(editInfor!=NULL)
	 size = editInfor->size();
	else size = 0;
	return size;
}



bool Alignment::insertSubstitution(int position, char base)
{
	if(this->editInfor==NULL)
		this->editInfor = new map<int, char>();
	this->editInfor->insert(std::pair<int, char>(position, base));
	return true;
}

SubjectEdge::SubjectEdge(SubjectRead* sRead)
{
	this->subjectRead = sRead;
	this->alignmentList = NULL;
	this->contained_alignmentList = NULL;
	this->DuplicateReadList = NULL;
	this->hashkeypositionmap = NULL;


}
SubjectEdge::~SubjectEdge()
{
	this->subjectRead = NULL;
	if(this->alignmentList!=NULL)
	{
		this->alignmentList->clear();
		delete this->alignmentList;
	}
	if(this->DuplicateReadList!=NULL)
	{
		this->DuplicateReadList->clear();
		delete this->DuplicateReadList;
	}
	if(this->hashkeypositionmap!=NULL)
	{
		for(std::map<UINT64, vector<int>*>::iterator it=hashkeypositionmap->begin(); it!=hashkeypositionmap->end(); ++it)
		{
			vector<int>* intvector = it->second;
			intvector->clear();
			delete intvector;
			it->second = NULL;
		}
		this->hashkeypositionmap->clear();
		delete this->hashkeypositionmap;
	}

}

void SubjectEdge::clearAlignmentList()
{
	if(this->alignmentList!=NULL)
	{
		this->alignmentList->clear();
		delete this->alignmentList;
		this->alignmentList = NULL;
	}

}

bool SubjectEdge::addKeyToMap(UINT64 keyvalue, int position)
{
	if(this->hashkeypositionmap==NULL)
		this->hashkeypositionmap = new map<UINT64, vector<int>*>();

	map<UINT64, vector<int>*>::const_iterator it = this->hashkeypositionmap->find(keyvalue);
	if(it!=this->hashkeypositionmap->end())
	{
		vector<int>* positionlist = it->second;
		positionlist->push_back(position);
	}
	else
	{
		vector<int>* positionlist = new vector<int>();
		positionlist->push_back(position);
		this->hashkeypositionmap->insert(std::pair<UINT64, vector<int>*>(keyvalue,positionlist));
	}

	return true;
}


bool SubjectEdge::addAlignment(Alignment* subjectAlignment)
{
	if(this->alignmentList==NULL)
		this->alignmentList = new vector<Alignment*>();
	this->alignmentList->push_back(subjectAlignment);
	return true;
}

bool SubjectEdge::addDuplicateList(QueryRead* queryRead)
{
	if(this->DuplicateReadList==NULL)
		this->DuplicateReadList = new vector<QueryRead*>();
	this->DuplicateReadList->push_back(queryRead);
	return true;
}
bool SubjectEdge::addContainedAlignment(ContainedAlignment* subjectAlignment)
{
	if(this->contained_alignmentList==NULL)
		this->contained_alignmentList = new vector<ContainedAlignment*>();
	this->contained_alignmentList->push_back(subjectAlignment);
	return true;
}

ContainedAlignment::ContainedAlignment(SubjectRead* sRead, QueryRead* qRead)
{
	// TODO Auto-generated constructor stub
	subjectRead = sRead;
	queryRead = qRead;

}

ContainedAlignment::~ContainedAlignment() {
	// TODO Auto-generated destructor stub

	queryRead = NULL;
	subjectRead = NULL;

}

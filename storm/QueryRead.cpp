/*
 * QueryRead.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "QueryRead.h"

QueryRead::QueryRead() {
	omp_init_lock(&lock);
	readSequence = "";
	readRevSequence = "";
	readName = "";
	readID = 0;
	frequency = 0;
	flag4Removal = false;
	queryAlignmentList=NULL;
	pairEndRead=NULL;
	pairPosition = 0;
	newalignmentcounter = 0;
}

QueryRead::QueryRead(string& sequence, string& name)
{
	omp_init_lock(&lock);
	readSequence = sequence;
	readRevSequence = "";
	readName = name;
	readID = 0;
	frequency = 1;
	flag4Removal = false;
	queryAlignmentList=NULL;
	pairEndRead=NULL;
	pairPosition = 0;
	newalignmentcounter = 0;
}


QueryRead::~QueryRead() {
	// TODO Auto-generated destructor stub
	if(queryAlignmentList!=NULL)
	{
	for(UINT64 i = 0; i< queryAlignmentList->size();i++)
		if(queryAlignmentList->at(i)!=NULL)
		delete queryAlignmentList->at(i);
	queryAlignmentList->clear();
	delete queryAlignmentList;
	queryAlignmentList=NULL;
	}

	//deletion of pairEndRead is in QueryDataset
	pairEndRead=NULL;
	omp_destroy_lock(&lock);

}

void QueryRead::setRemovalTag(bool tag)
{
//	omp_set_lock(&lock);
#pragma omp atomic
	flag4Removal &= tag;
//	omp_unset_lock(&lock);
}

bool QueryRead::addAlignmentWithLock(Alignment* subjectAlignment)
{
	omp_set_lock(&lock);
	if(this->queryAlignmentList==NULL)
	{
		this->queryAlignmentList = new vector<Alignment*>();
//		this->queryAlignmentList->reserve(MAX_VECTOR);
	}
	this->queryAlignmentList->push_back(subjectAlignment);
	this->newalignmentcounter++;
	omp_unset_lock(&lock);
	return true;
}

bool QueryRead::needRemoval()
{
	return flag4Removal;
}
void QueryRead::setName(string & name)
{
	readName = name;
}
void QueryRead::setSequence(string & sequence)
{
	readSequence = sequence;
}

void QueryRead::setSequence(string & sequence, string & reverseSequence)
{
	readSequence = sequence;
	readRevSequence = reverseSequence;

}

void QueryRead::setFrequency(UINT32 number)
{
	frequency = number;
}
void QueryRead::setIdentifier(UINT64 id)
{
	readID = id;
}
string QueryRead::getName()
{
	return readName;
}
string QueryRead::getSequence()
{
	return readSequence;
}

UINT64 QueryRead::getIdentifier()
{
	return readID;
}
UINT32 QueryRead::getFrequency()
{
	return frequency;
}
UINT32 QueryRead::getReadLength()
{
	return this->readSequence.length();
}

UINT64 QueryRead::getAlignmentNumber()
{
	if(this->queryAlignmentList!=NULL)
		return this->queryAlignmentList->size();
	else
		return 0;
}

bool QueryRead::addAlignment(Alignment* subjectAlignment)
{
	if(this->queryAlignmentList==NULL)
	{
		this->queryAlignmentList = new vector<Alignment*>();
//		this->queryAlignmentList->reserve(MAX_VECTOR);
	}
	this->queryAlignmentList->push_back(subjectAlignment);
	this->newalignmentcounter++;
	return true;
}
string QueryRead::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C <==> G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}
string QueryRead::reverseComplement()
{
	if(this->readRevSequence=="")
		this->readRevSequence = this->QueryRead::reverseComplement(this->readSequence);
	return this->readRevSequence;
}

bool wayToSort(int i, int j) { return i < j; }

bool QueryRead::printAlignmentToFile(ofstream & filePointer)
{

	if(this->queryAlignmentList==NULL)return false;
	for(UINT64 i =0; i< this->queryAlignmentList->size();i++)
	{

	Alignment *align = this->queryAlignmentList->at(i);
	if(align!=NULL)
	{
	string outputString = this->alignmentToOmegaString(align);
//	cout<<outputString;
	filePointer<<outputString;
	}
	}

	return true;
}


string QueryRead::getStringForPrintAlignmentToFile()
{
	string outputString = "";

	if(this->queryAlignmentList==NULL)return outputString;
	std::stringstream sstm;
	sstm<<"";
	for(UINT64 i =0; i< this->queryAlignmentList->size();i++)
	{

		Alignment *align = this->queryAlignmentList->at(i);
		if(align!=NULL)
		{
		string currentString = this->alignmentToOmegaString(align);
		sstm<<currentString;
		}

	}

	outputString = sstm.str();
//	cout<<outputString;

	return outputString;
}

string QueryRead::getStringForPrintAlignmentToFileWithoutDuplicate()
{
	string outputString = "";

	if(this->queryAlignmentList==NULL)return outputString;
	std::stringstream sstm;
	sstm<<"";
	for(UINT64 i =0; i< this->queryAlignmentList->size();i++)
	{

	Alignment *align = this->queryAlignmentList->at(i);
	if(align!=NULL && align->subjectRead->getName()>align->queryRead->getName())
	{
	string currentString = this->alignmentToOmegaString(align);
	sstm<<currentString;
	}

	}

	outputString = sstm.str();
//	cout<<outputString;

	return outputString;
}

string QueryRead::alignmentToOmegaString(Alignment* align)
{
	string outputString="";
	// query print first for good orientation
	int orient = align->orientationTranslate();

	if(orient<0) return outputString;
	/*SourceVertexId       DestinationVertexId     Properties

	where properties should have

	orientation, overlap length, substitutions, edits, length1, start1,
	stop1, length2, start2, stop2, error info
	*/
	vector<int> coordinates;
	coordinates.push_back(align->queryEnd);
	coordinates.push_back(align->subjectStart);
	coordinates.push_back(align->subjectEnd);
	coordinates.push_back(0);
	sort(coordinates.begin(), coordinates.end(), wayToSort);
	int overlapLength = coordinates.at(2)-coordinates.at(1)+1;
	int queryStart = coordinates.at(1);
	int queryEnd = coordinates.at(2);
//	int subStart = coordinates.at(1)-align->subjectStart;
//	int subEnd = coordinates.at(2)-align->subjectStart;
	int subjectOrigin = align->subjectStart;
	std::stringstream sstm;

	string errorinfo = align->editInfoToOmegaString();
	int editdistance = align->getEditDistance();
	int numberofsubsitution = editdistance;

	//change coordinates if the alignment stretching to the left side
	if(align->subjectStart<0)
	{
		int tempqueryStart=queryStart;
		queryStart = align->queryEnd-queryEnd;
		queryEnd = align->queryEnd-tempqueryStart;
		subjectOrigin = align->queryEnd - align->subjectEnd;
	}

	int subStart = queryStart-subjectOrigin;
	int subEnd = queryEnd - subjectOrigin;
	sstm <<  this->readName << "\t" << align->subjectRead->getName() << "\t" <<orient <<
			"," << overlapLength << "," << numberofsubsitution << "," << editdistance
			<< "," << this->getReadLength() << "," << queryStart << "," << queryEnd <<","
			<< align->subjectRead->getReadLength() << "," << subStart<< "," << subEnd<<","<<errorinfo<<endl;

	outputString = sstm.str();

	return outputString;

}

PairedEndRead* QueryRead::getPairEndRead()
{
	return this->pairEndRead;
}
void QueryRead::setPairEndRead(PairedEndRead* pairedEndRead)
{
	this->pairEndRead = pairedEndRead;
}
UINT8 QueryRead::getPairPosition()
{
	return this->pairPosition;
}
void QueryRead::setPairPosition(UINT8 pos)
{
	this->pairPosition = pos;
}

bool QueryRead::addRightMatePair(QueryRead *rightRead, UINT8 orientation)
{
	if(rightRead->getPairPosition()==0)
	{
		PairedEndRead * pair = new PairedEndRead();
		pair->leftRead = this;
		pair->rightRead = rightRead;
		pair->orientation = orientation;
		this->setPairEndRead(pair);
		this->setPairPosition(1);
		rightRead->setPairEndRead(pair);
		rightRead->setPairPosition(2);
		return true;
	}
	else if(rightRead->getPairPosition()==1)return false;
	else if(rightRead->getPairPosition()==2 && rightRead->getPairEndRead()->leftRead->getName()!=this->getName())return false;
	else return true;

}

bool compareRightAlignments (Alignment *align1, Alignment *align2)
{
	int align1_start, align2_start;
	if(align1->queryOrientation==true)
		align1_start = align1->subjectStart;
	else
		align1_start = align1->queryEnd - align1->subjectEnd;

	if(align2->queryOrientation==true)
		align2_start = align2->subjectStart;
	else
		align2_start = align2->queryEnd - align2->subjectEnd;

	return align1_start < align2_start;
}
void QueryRead::sortRightAlignments(vector<Alignment*>* alignList)
{
	if(alignList!=NULL)
	std::sort(alignList->begin(),alignList->end(), compareRightAlignments);
	else
	cout<<"no reads in the dataset";

}

bool compareLeftAlignments (Alignment *align1, Alignment *align2)
{
	int align1_start, align2_start;
	if(align1->queryOrientation==false)
		align1_start = align1->subjectStart;
	else
		align1_start = align1->queryEnd - align1->subjectEnd;

	if(align2->queryOrientation==false)
		align2_start = align2->subjectStart;
	else
		align2_start = align2->queryEnd - align2->subjectEnd;

	return align1_start < align2_start;
}
void QueryRead::sortLeftAlignments(vector<Alignment*>* alignList)
{
	if(alignList!=NULL)
	std::sort(alignList->begin(),alignList->end(), compareLeftAlignments);
	else
	cout<<"no reads in the dataset";

}

string QueryRead::getStringForPrintNonTransitiveAlignmentToFile()
{
	string outputString = "";
	if(this->queryAlignmentList==NULL||this->queryAlignmentList->size()==0)return outputString;


	this->markTransitiveEdge();

	std::stringstream sstm;
	sstm<<"";
	for(UINT64 i =0; i< this->queryAlignmentList->size();i++)
	{
		Alignment *align = this->queryAlignmentList->at(i);
		if(align!=NULL&&align->containedTag!=128&&align->subjectRead->getName()>align->queryRead->getName())
		{

			string alignstring = this->alignmentToOmegaString(align);
			sstm<<alignstring;
		}
	}
	outputString = sstm.str();
	return outputString;
}



bool QueryRead::markTransitiveEdgeIncremental()
{
	if(this->queryAlignmentList==NULL||this->queryAlignmentList->size()==0) return false;
	INT64 lastIndex = this->queryAlignmentList->size()-1;
	INT64 firstNewIndex = lastIndex - this->newalignmentcounter +1;

	for(INT64 i=firstNewIndex; i<=lastIndex; i++)
	{
		Alignment * newAlign = this->queryAlignmentList->at(i);
		if(newAlign->containedTag!=0)continue;
		bool isNewAlignRight; int newAlignStart;
		newAlign->isLeftOrRightAlignment(isNewAlignRight, newAlignStart);
		for(INT64 j=i-1; j>=0; j--)
		{
			Alignment * oldAlign = this->queryAlignmentList->at(j);
			bool isOldAlignRight; int oldAlignStart;
			oldAlign->isLeftOrRightAlignment(isOldAlignRight, oldAlignStart);
			if(isNewAlignRight==isOldAlignRight)
			{
				Alignment* inner; Alignment *outter;
				if(newAlignStart>oldAlignStart)
				{
					inner = oldAlign;
					outter = newAlign;
				}
				else if(newAlignStart<oldAlignStart)
				{
					inner = newAlign;
					outter = oldAlign;
				}
				else if(newAlign->subjectRead->getReadLength()<oldAlign->subjectRead->getReadLength())
				{
					inner = oldAlign;
					outter = newAlign;
				}
				else
				{
					inner = newAlign;
					outter = oldAlign;
				}
				//there is no need to comparison if outter is already marked
				if(outter->containedTag!=128)
				{
					bool alignsuccess;
					if(isNewAlignRight==true)
						alignsuccess = this->checkRightAlignmentOverlap(outter, inner);
					else
						alignsuccess = this->checkLeftAlignmentOverlap(outter, inner);
					if(alignsuccess)
					{
						outter->containedTag = 128;
						outter->subjectRead->decreaselinknumberWithLock();
					}
				}
			}
		}
	}

	//reset the counter after using it
	this->newalignmentcounter = 0;
	return true;
}

bool QueryRead::markTransitiveEdge()
{
	vector<Alignment* > * rightAlignList = new vector<Alignment* >();
	vector<Alignment* > * leftAlignList = new vector<Alignment* >();

	if(this->queryAlignmentList==NULL)return false;
	for(UINT64 i=0; i<this->queryAlignmentList->size(); i++)
	{
		Alignment* align = this->queryAlignmentList->at(i);
		if(align->containedTag==0)
		{
		if(align->queryOrientation==true)
		{
			if(align->subjectStart>=0)
				rightAlignList->push_back(align);
			else
				leftAlignList->push_back(align);
		}
		else
		{
			if(align->subjectStart>=0)
				leftAlignList->push_back(align);
			else
				rightAlignList->push_back(align);
		}
		}
	}



	if(rightAlignList->size()>1)
	{
		this->sortRightAlignments(rightAlignList);
		for(int i=rightAlignList->size()-1;i>0;i--)
		{
			Alignment* targetAlign = rightAlignList->at(i);
			for(int j=i-1;j>=0;j--)
			{
				Alignment* interiorAlign = rightAlignList->at(j);
				if(this->checkRightAlignmentOverlap(targetAlign,interiorAlign))
				{
					targetAlign->containedTag = 128;
					targetAlign->subjectRead->decreaselinknumberWithLock();
					break;
				}
			}
		}

	}

	if(leftAlignList->size()>1)
	{
		this->sortLeftAlignments(leftAlignList);
		for(int i=leftAlignList->size()-1;i>0;i--)
		{
			Alignment* targetAlign = leftAlignList->at(i);
			for(int j=i-1;j>=0;j--)
			{
				Alignment* interiorAlign = leftAlignList->at(j);
				if(this->checkLeftAlignmentOverlap(targetAlign,interiorAlign))
				{
					targetAlign->containedTag = 128;
					targetAlign->subjectRead->decreaselinknumberWithLock();
					break;
				}
			}
		}
	}

	rightAlignList->clear();
	delete rightAlignList;
	leftAlignList->clear();
	delete leftAlignList;
	return true;
}

void QueryRead::deleteTransitiveRemovableRedundantEdge()
{

if(this->queryAlignmentList!=NULL&&this->queryAlignmentList->size()>0)
{
	INT64 topElement = 0;
	INT64 bottomElement = (int)this->queryAlignmentList->size()-1;
	UINT64 nonEmptySpace = this->queryAlignmentList->size();
	//bottomElement needs to point to the first good element from the bottom
	while(bottomElement>=0&&this->queryAlignmentList->at(bottomElement)->containedTag>0)
	{
		delete this->queryAlignmentList->at(bottomElement);
		nonEmptySpace = nonEmptySpace -1 ;
		bottomElement--;
	}
	while(bottomElement>topElement)
	{
		//check top element
		if(this->queryAlignmentList->at(topElement)->containedTag>0)
		{
			delete this->queryAlignmentList->at(topElement);
			nonEmptySpace = nonEmptySpace -1 ;
			this->queryAlignmentList->at(topElement) = this->queryAlignmentList->at(bottomElement);
			bottomElement--;
			topElement++;
			while(bottomElement>=topElement&&this->queryAlignmentList->at(bottomElement)->containedTag>0)
			{
				delete this->queryAlignmentList->at(bottomElement);
				bottomElement--;
				nonEmptySpace = nonEmptySpace -1 ;
			}
		}
		else
		{
			topElement++;
		}
	}

	this->queryAlignmentList->resize(nonEmptySpace);

}

}

bool QueryRead::checkRightAlignmentOverlap(Alignment* targetAlign,Alignment* interiorAlign)
{

	if(targetAlign->queryOrientation==false && interiorAlign->queryOrientation == false)
	{
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = targetAlign->subjectEnd;
		start = (targetAlign->subjectStart>=interiorAlign->subjectStart)?targetAlign->subjectStart:interiorAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(start,-1);
				string interiorsequence = interiorAlign->getSubjectSubSequence(start,-1);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else if(targetAlign->queryOrientation==true && interiorAlign->queryOrientation==true)
	{
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = (targetAlign->subjectEnd<=interiorAlign->subjectEnd)?targetAlign->subjectEnd:interiorAlign->subjectEnd;
		start = targetAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				string interiorsequence = interiorAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else if(targetAlign->queryOrientation==true && interiorAlign->queryOrientation==false)
	{
		interiorAlign->reverseAlignment();
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = (targetAlign->subjectEnd<=interiorAlign->subjectEnd)?targetAlign->subjectEnd:interiorAlign->subjectEnd;
		start = targetAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				string interiorsequence = interiorAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else if(targetAlign->queryOrientation==false && interiorAlign->queryOrientation==true)
	{
		targetAlign->reverseAlignment();
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = (targetAlign->subjectEnd<=interiorAlign->subjectEnd)?targetAlign->subjectEnd:interiorAlign->subjectEnd;
		start = targetAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				string interiorsequence = interiorAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else return false;
}


bool QueryRead::checkLeftAlignmentOverlap(Alignment* targetAlign,Alignment* interiorAlign)
{

	if(targetAlign->queryOrientation==true && interiorAlign->queryOrientation == true)
	{
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = targetAlign->subjectEnd;
		start = (targetAlign->subjectStart>=interiorAlign->subjectStart)?targetAlign->subjectStart:interiorAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(start,-1);
				string interiorsequence = interiorAlign->getSubjectSubSequence(start,-1);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else if(targetAlign->queryOrientation==false && interiorAlign->queryOrientation==false)
	{
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = (targetAlign->subjectEnd<=interiorAlign->subjectEnd)?targetAlign->subjectEnd:interiorAlign->subjectEnd;
		start = targetAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				string interiorsequence = interiorAlign->getSubjectSubSequence(targetAlign->queryEnd+1,end);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else if(targetAlign->queryOrientation==false && interiorAlign->queryOrientation==true)
	{
		targetAlign->reverseAlignment();
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = targetAlign->subjectEnd;
		start = (targetAlign->subjectStart>=interiorAlign->subjectStart)?targetAlign->subjectStart:interiorAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(start,-1);
				string interiorsequence = interiorAlign->getSubjectSubSequence(start,-1);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else if(targetAlign->queryOrientation==true && interiorAlign->queryOrientation==false)
	{
		interiorAlign->reverseAlignment();
		int totalErrors=0;
		int start, end, overlaplength;
		int maxError;
		end = targetAlign->subjectEnd;
		start = (targetAlign->subjectStart>=interiorAlign->subjectStart)?targetAlign->subjectStart:interiorAlign->subjectStart;
		overlaplength = end - start +1;
		if(overlaplength>=Config::minimumOverlapLength)
		{
			if(Config::useErrorRate==true)
			{
				maxError = floor(overlaplength*((double)(Config::maxErrorRate))/100);

			}
			else maxError = Config::maxMismatch;
			totalErrors = targetAlign->distanceBetweenEditInfo(interiorAlign);
			if(totalErrors>maxError)return false;
			else
			{
				string targetsequence = targetAlign->getSubjectSubSequence(start,-1);
				string interiorsequence = interiorAlign->getSubjectSubSequence(start,-1);
				for(int i =0; i<(int)targetsequence.length();i++)
				{
					if(targetsequence.at(i)!=interiorsequence.at(i))
						totalErrors++;
					if(totalErrors>maxError)return false;
				}
				return true;
			}
		}
		else return false;
	}
	else return false;
}
/*
bool QueryRead::correctErrors(int minimumDepth, double cutoff){

	for(int i = 0; i < queryAlignmentList.size(); i++){

	}
}
*/


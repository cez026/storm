/*
 * PairEndReadMerger.cpp
 *
 *  Created on: Jun 29, 2015
 *      Author: qy2
 */

#include "PairEndReadMerger.h"
#define MINIMUMDEPTH 5
#define RECALLPERCENT 0.9

PairEndReadMerger::PairEndReadMerger(HashTableMethod * hashTableMethod) {


	this->hashTableMethod = hashTableMethod;

	this->left_minimumOverlapLength = Config::left_minimumOverlapLength;
	this->right_minimumOverlapLength = Config::right_minimumOverlapLength;
	this->overall_minimumOverlapLength = Config::overall_minimumOverlapLength;

	if(Config::hashtabletype=="single")
		this->checkKeyLength = Config::hashKeyLength;
	else
		this->checkKeyLength = Config::hashKeyLength_left;


	this->maxErrorRate = ((double)(Config::maxErrorRate))/100;
	this->maxMismatch = Config::maxMismatch;
	this->minimumdepth = Config::minimumdepth;
	this->recallpercent = Config::recallpercent;

	if(Config::useErrorRate==true)
	{
		if(this->maxErrorRate==0)
			this->allowmismatch = false;
		else
			this->allowmismatch = true;
	}
	else
	{
		if(this->maxMismatch==0)
			this->allowmismatch = false;
		else
			this->allowmismatch = true;
	}

}

PairEndReadMerger::~PairEndReadMerger() {

}

bool PairEndReadMerger::checkRightAlignment(string subjectSequence,string rightQuerySequence,int alignIndex, bool direction, int tolerateError)
{
	int currentError = 0;
	string right,subject;
	bool successflag = true;
	if(direction==true)
	{
		int stringlen = ((int)rightQuerySequence.length()<(int)subjectSequence.length()-alignIndex)?(int)rightQuerySequence.length():(int)subjectSequence.length()-alignIndex;
//		if(stringlen-1>=(int)rightQuerySequence.length())cout<<"1: substr(): "<<endl;//SUBSTR
		right = rightQuerySequence.substr(0,stringlen);
//		if(alignIndex<0||alignIndex+stringlen-1>=(int)subjectSequence.length())cout<<"2: substr(): "<<endl;//SUBSTR
		subject = subjectSequence.substr(alignIndex,stringlen);
	}
	else
	{
		int stringlen = ((int)rightQuerySequence.length()<1+alignIndex)?(int)rightQuerySequence.length():1+alignIndex;
//		if((int)rightQuerySequence.length()-stringlen<0||(int)rightQuerySequence.length()-stringlen, stringlen-1>=(int)rightQuerySequence.length())cout<<"3: substr(): "<<endl;//SUBSTR
		right = rightQuerySequence.substr((int)rightQuerySequence.length()-stringlen, stringlen);
//		if(alignIndex+1-stringlen<0||alignIndex+1-stringlen+stringlen-1>=(int)subjectSequence.length())cout<<"4: substr(): "<<endl;//SUBSTR
		subject = subjectSequence.substr(alignIndex+1-stringlen,stringlen);
	}
	int cursor = 0;
	while(successflag==true&&cursor<(int)right.length())
	{
		if(currentError==tolerateError)
		{
//			successflag = (right.substr(cursor,right.length()-cursor)==subject.substr(cursor,subject.length()-cursor));
			successflag = (strcmp(right.substr(cursor,right.length()-cursor).c_str(),subject.substr(cursor,subject.length()-cursor).c_str())==0);
			cursor=right.length();
		}
		else if(currentError<tolerateError)
		{
			if(right.at(cursor)!=subject.at(cursor))
				currentError++;
			cursor++;
		}
		else
		{
			successflag = false;
		}
	}
	return successflag;
}

bool PairEndReadMerger::createPairedEndAlignment(map<UINT64, vector<int>*>* hashkeypositionmap, Alignment* alignment, QueryRead* rightqueryread,UINT16 keylength, UINT16 rightminimumoverlap, UINT16 overallminimumoverlap)
{

	if(rightqueryread->getSequence().length()<rightminimumoverlap||alignment->subjectRead->getSequence().length()<rightminimumoverlap) return false;
	bool direction = alignment->queryOrientation;
	vector<int>* positionlist=NULL;
	int offset = 0;//offset for the start point (direction = true) or the end point (direction = false)
	string rightquerysequence = "";
	if(direction==true)
	{
		//check alignment with query 2, with query 1 hit by hash key
		//query 1: XXXXXXXXXXXKKKKKK query 2: MMMMMMMMXXXXXXXXX
		//subject:            KKKKKKXXXXXXXXXXMMMMMMMM
		rightquerysequence = rightqueryread->getSequence();
//		if(keylength-1>=(int)rightquerysequence.length())cout<<"5: substr(): "<<endl;//SUBSTR
		string rightkeystring = rightquerysequence.substr(0,keylength);
		UINT64 rightkeyvalue = this->hashTableMethod->getHashTable()->hashFunction(rightkeystring);

		map<UINT64, vector<int>*>::const_iterator it = hashkeypositionmap->find(rightkeyvalue);
		if(it!=hashkeypositionmap->end())
		{
			positionlist = it->second;
			offset = 0;
		}
		else if(this->allowmismatch)
		{
//			if(rightminimumoverlap-keylength<0||rightminimumoverlap-keylength+keylength-1>=(int)rightquerysequence.length())cout<<"6: substr(): "<<endl;//SUBSTR
			string alternativerightkeystring = rightquerysequence.substr(rightminimumoverlap-keylength, keylength);
			UINT64 alternativerightkeyvalue = this->hashTableMethod->getHashTable()->hashFunction(alternativerightkeystring);
			map<UINT64, vector<int>*>::const_iterator it2 = hashkeypositionmap->find(alternativerightkeyvalue);
			if(it2!=hashkeypositionmap->end())
			{
				positionlist = it2->second;
				offset= -(rightminimumoverlap-keylength);
			}
		}

	}
	else
	{
		//check alignment with query 2 (reverse), with query 1 (reverse) hit by hash key
		//query 2 (rev): XXXXXXXXXXXMMMMMM query 1 (rev): KKKKKKKKXXXXXXXXX
		//subject (fwd):            MMMMMMXXXXXXXXXXXXXXXXKKKKKKKK
		rightquerysequence = rightqueryread->reverseComplement();
//		if((int)rightquerysequence.length()-keylength<0||(int)rightquerysequence.length()-keylength+keylength-1>=(int)rightquerysequence.length())cout<<"7: substr(): "<<endl;//SUBSTR
		string rightkeystring = rightquerysequence.substr(rightquerysequence.length()-keylength,keylength);
		UINT64 rightkeyvalue = this->hashTableMethod->getHashTable()->hashFunction(rightkeystring);

		map<UINT64, vector<int>*>::const_iterator it = hashkeypositionmap->find(rightkeyvalue);
		if(it!=hashkeypositionmap->end())
		{
			positionlist = it->second;
			offset = keylength-1;
		}
		else if(this->allowmismatch)
		{
//			if((int)rightquerysequence.length()-rightminimumoverlap<0||(int)rightquerysequence.length()-rightminimumoverlap+keylength-1>=(int)rightquerysequence.length())cout<<"8: substr(): "<<endl;//SUBSTR
			string alternativerightkeystring = rightquerysequence.substr(rightquerysequence.length()-rightminimumoverlap, keylength);
			UINT64 alternativerightkeyvalue = this->hashTableMethod->getHashTable()->hashFunction(alternativerightkeystring);
			map<UINT64, vector<int>*>::const_iterator it2 = hashkeypositionmap->find(alternativerightkeyvalue);
			if(it2!=hashkeypositionmap->end())
			{
				positionlist = it2->second;
				offset= rightminimumoverlap-1;
			}
		}
	}

	if(positionlist==NULL) return false;

	bool successflag = false;
	for(int i=0; i< (int)positionlist->size(); i++)
	{
		int pairedEndAlign = positionlist->at(i);
		pairedEndAlign = pairedEndAlign+offset;
		string subjectSequence = alignment->subjectRead->getSequence();
		if(direction==true&&(pairedEndAlign<0||pairedEndAlign+rightminimumoverlap>(int)subjectSequence.length()))
			continue;
		if(direction==false&&(pairedEndAlign<rightminimumoverlap-1||pairedEndAlign>=(int)subjectSequence.length()))
			continue;

		 int leftqueryend,subjectstart,subjectend, rightquerystart, rightqueryend;
		 bool suborient;
		 if(alignment->queryOrientation==true)
		 {
			 leftqueryend = alignment->queryEnd;
			 subjectstart = alignment->subjectStart;
			 subjectend = alignment->subjectEnd;
			 rightquerystart = pairedEndAlign+subjectstart;
			 rightqueryend = rightquerystart-1+rightqueryread->getReadLength();
			 suborient = true;

		 }
		 else
		 {
			 leftqueryend = alignment->queryEnd;
			 subjectstart = -(alignment->subjectEnd)+alignment->queryEnd;
			 subjectend = -(alignment->subjectStart)+alignment->queryEnd;
			 rightquerystart = -(pairedEndAlign+alignment->subjectStart) + alignment->queryEnd;
			 rightqueryend = rightquerystart-1+rightqueryread->getReadLength();
			 suborient = false;

		 }
		 if(rightquerystart>0&&rightqueryend>leftqueryend&&leftqueryend-subjectstart+1+subjectend-rightquerystart+1>=overallminimumoverlap)
		 {
			 int lefterror = alignment->getEditDistance();
			 int totalerror=0;
			 if(this->allowmismatch)
			 {
				if(Config::useErrorRate==true)
				{
					totalerror = floor((subjectend-rightquerystart+1)*this->maxErrorRate);
				}
				else
				{
					totalerror = this->maxMismatch;

				}
			 }
			 int righterror=totalerror-lefterror;
			 if(checkRightAlignment(alignment->subjectRead->getSequence(),rightquerysequence,pairedEndAlign, direction, righterror))
			 {
				string subjectSubString;
				string pre;
				int interval = rightquerystart-leftqueryend-1;
				int start, end;
				if(interval<0)
				{
					start = rightquerystart;
					end = leftqueryend;
					pre="-";
				}
				else
				{
					start = leftqueryend+1;
					end = rightquerystart-1;
					pre="+";
				}

				if(suborient==true)
				{
//					if(start-subjectstart<0||start-subjectstart+end-start+1-1>=(int)alignment->subjectRead->getSequence().length())cout<<"9: substr(): "<<endl;//SUBSTR
					subjectSubString = alignment->subjectRead->getSequence().substr(start-subjectstart, end-start+1);

				}
				else
				{
//					if(subjectend-end<0||subjectend-end+end-start+1-1>=(int)alignment->subjectRead->getSequence().length())cout<<"10: substr(): "<<endl;//SUBSTR
					string temp = alignment->subjectRead->getSequence().substr(subjectend-end, end-start+1);
					subjectSubString = alignment->subjectRead->reverseComplement(temp);
				}
				string gapstring = pre+subjectSubString;

//				#pragma omp critical(addGapStringToPairEndReads)
//				 {
				 alignment->queryRead->getPairEndRead()->insertgapstringfrequencymap(gapstring);
//				 }
//				 successflag = true;
			 }
		 }

	}//end of for
	return successflag;

}

/*
bool PairedEndReadsMerger::start()
{
	if(Config::useErrorRate==true)
		return this->start_errorrate_seqfree();
	else
		return this->start_errornum_seqfree();
	return false;
}
*/

bool PairEndReadMerger::start()
{

	CLOCKSTART;
	MEMORYSTART;


	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	if(this->left_minimumOverlapLength>this->right_minimumOverlapLength)
		subjectDataset->setReadFilterSize(this->left_minimumOverlapLength);
	else
		subjectDataset->setReadFilterSize(this->right_minimumOverlapLength);

	vector<SubjectEdge*>* subjectEdgesList = new  vector<SubjectEdge*>();
	subjectEdgesList->clear();
	while(subjectDataset->loadNextChunk(subjectEdgesList))
	{

omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgesList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgesList->at(i);
			this->hashTableMethod->searchHashTableFullScanAndRecordKeys(subjectEdge); //need a full scan here

			if(subjectEdge->alignmentList!=NULL)
			{
				map<UINT64, vector<int>*>* hashkeypositionmap = subjectEdge->hashkeypositionmap;
//				cout<<subjectEdge->subjectRead->getName()<<"  "<<subjectEdge->alignmentList->size()<<endl;
				//add alignments to the query reads before the edges are destroyed.
				for(INT64 j = 0; j < (int)subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);

					 QueryRead * rightQueryRead = alignment->queryRead->getPairEndRead()->rightRead;
					 this->createPairedEndAlignment(hashkeypositionmap, alignment,rightQueryRead,this->checkKeyLength, this->right_minimumOverlapLength, this->overall_minimumOverlapLength);

					 delete alignment;//this alignment need to be deleted.
				}

			}
			delete subjectEdge->subjectRead;
			delete subjectEdge;
			subjectEdgesList->at(i)=NULL;
		}

	}	// end of omp parallel

		//clear the subject read vector for the next round
		subjectEdgesList->clear();
		subjectEdgesList->resize(0);

	}// end of while chunk loop

	{//BLOCK of merging
		cout<<"this is merging process----"<<endl;
		CLOCKSTART;
		MEMORYSTART;
	//need to call the base pairs per end reads pair
	vector<PairedEndRead*>* pairedReadList = this->hashTableMethod->getDataset()->getPairedEndReadList();

omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(int i=0; i<(int)pairedReadList->size();i++)
		{
			PairedEndRead* pairedEndRead = pairedReadList->at(i);
			pairedEndRead->reconstructSequencefromGapSeq(this->minimumdepth, this->recallpercent);
		}
	}

	MEMORYSTOP;
	CLOCKSTOP;
	}//end of BLOCK of merging

	string notMergedFileName = Config::outputfilename+"notmerged.fasta";
	string mergedFileName = Config::outputfilename+"merged.fasta";
	printMergedToFile(notMergedFileName, mergedFileName);

	delete subjectEdgesList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;

	return true;

}

bool PairEndReadMerger::printMergedToFile(string notMergedFileName, string mergedFileName)
{
	CLOCKSTART;
	MEMORYSTART;
	ofstream mergeFilePointer;
	ofstream notMergeFilePointer;
	mergeFilePointer.open(mergedFileName.c_str());
	notMergeFilePointer.open(notMergedFileName.c_str());
	if(mergeFilePointer == NULL || notMergeFilePointer ==NULL)
	{
		cout<<"Unable to open file: "<<endl;
		return false;
	}
	vector<PairedEndRead*>* pairedReadList = this->hashTableMethod->getDataset()->getPairedEndReadList();


	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
	#pragma omp parallel shared(mergeFilePointer,notMergeFilePointer)
		{
		#pragma omp for schedule(dynamic)
			for(int i=0; i<(int)pairedReadList->size();i++)
			{
				PairedEndRead* pairedEndRead = pairedReadList->at(i);
				if(pairedEndRead->mergedSequence.length()>0)
				{
					std::stringstream sstm;
					sstm<<">"<<"LENGTH="<<pairedEndRead->mergedSequence.length()<<"&&&&"<<pairedEndRead->leftRead->getName()<<"&&&&"<<pairedEndRead->rightRead->getName()<<endl;
					sstm<<pairedEndRead->mergedSequence<<endl;
					string output = sstm.str();
		#pragma omp critical(printToMerged)
					mergeFilePointer<<output;
				}
				else
				{
					std::stringstream sstm;
					switch(Config::pairedEndReadsInputMode)
					{
					case 0:
						sstm<<">"<<pairedEndRead->leftRead->getName()<<endl;
						sstm<<pairedEndRead->leftRead->reverseComplement()<<endl;
						sstm<<">"<<pairedEndRead->rightRead->getName()<<endl;
						sstm<<pairedEndRead->rightRead->reverseComplement()<<endl;
						break;
					case 1:
						sstm<<">"<<pairedEndRead->leftRead->getName()<<endl;
						sstm<<pairedEndRead->leftRead->reverseComplement()<<endl;
						sstm<<">"<<pairedEndRead->rightRead->getName()<<endl;
						sstm<<pairedEndRead->rightRead->getSequence()<<endl;
						break;
					case 2:
						sstm<<">"<<pairedEndRead->leftRead->getName()<<endl;
						sstm<<pairedEndRead->leftRead->getSequence()<<endl;
						sstm<<">"<<pairedEndRead->rightRead->getName()<<endl;
						sstm<<pairedEndRead->rightRead->reverseComplement()<<endl;
						break;
					case 3:
						sstm<<">"<<pairedEndRead->leftRead->getName()<<endl;
						sstm<<pairedEndRead->leftRead->getSequence()<<endl;
						sstm<<">"<<pairedEndRead->rightRead->getName()<<endl;
						sstm<<pairedEndRead->rightRead->getSequence()<<endl;
						break;
					default:
						sstm<<">"<<pairedEndRead->leftRead->getName()<<endl;
						sstm<<pairedEndRead->leftRead->getSequence()<<endl;
						sstm<<">"<<pairedEndRead->rightRead->getName()<<endl;
						sstm<<pairedEndRead->rightRead->reverseComplement()<<endl;
					}

					string output = sstm.str();
		#pragma omp critical(printToNOTMerged)
					notMergeFilePointer<<output;
				}

			}

		}




	mergeFilePointer.close();
	notMergeFilePointer.close();
	MEMORYSTOP;
	CLOCKSTOP;
	return true;
}

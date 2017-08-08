/*
 * PairedEndReadsMerger.cpp
 *
 *  Created on: Mar 29, 2015
 *      Author: qy2
 */

#include "PairedEndReadsMerger.h"
#define MINIMUMDEPTH 5
#define RECALLPERCENT 0.9

PairedEndReadsMerger::PairedEndReadsMerger(HashTableMethod * hashTableMethod) {
	// TODO Auto-generated constructor stub

	this->hashTableMethod = hashTableMethod;
	this->left_minimumOverlapLength = Config::left_minimumOverlapLength;
	this->right_minimumOverlapLength = Config::right_minimumOverlapLength;
	this->overall_minimumOverlapLength = Config::overall_minimumOverlapLength;

	this->maxErrorRate = ((double)(Config::maxErrorRate))/100;
	this->maxMismatch = Config::maxMismatch;
}

PairedEndReadsMerger::~PairedEndReadsMerger() {
	// TODO Auto-generated destructor stub
}

//get the offest referred to the subject (forward) start point.
int PairedEndReadsMerger::checkAlignment(string subjectString, string queryString, UINT16 rightminimumoverlap, UINT16 overallminimumoverlap,int tolerateErrors, bool direction)
{

	if(direction == true)
	{
		//check alignment with query 2, with query 1 hit by hash key
		//query 1: XXXXXXXXXXXKKKKKK query 2: MMMMMMMMXXXXXXXXX
		//subject:            KKKKKKXXXXXXXXXXMMMMMMMM
		if(queryString.length()<rightminimumoverlap||subjectString.length()<rightminimumoverlap) return -1;
		for(int i=subjectString.length()-rightminimumoverlap;i>=0;i--)
		{
			int errorCount = 0;
			for(int j=i;j<(int)subjectString.length() && j-i<(int)queryString.length() && errorCount<=tolerateErrors;j++)
			{

				if(subjectString.at(j)!=queryString.at(j-i))errorCount++;

			}
			if(errorCount<=tolerateErrors) return i;
		}
		return -1;
	}
	else
	{
		//check alignment with query 2 (reverse), with query 1 (reverse) hit by hash key
		//query 2 (rev): XXXXXXXXXXXMMMMMM query 1 (rev): KKKKKKKKXXXXXXXXX
		//subject (rev):            MMMMMMXXXXXXXXXXXXXXXXKKKKKKKK
		if(queryString.length()<rightminimumoverlap||subjectString.length()<rightminimumoverlap) return -1;
		for(int i=rightminimumoverlap-1;i<(int)subjectString.length();i++)
		{
			int errorCount = 0;
			for(int j=i;j>=0&&i-j<(int)queryString.length()&& errorCount<=tolerateErrors;j--)
			{

				if(subjectString.at(j)!=queryString.at(queryString.length()-i-1+j))errorCount++;

			}
			if(errorCount<=tolerateErrors) return i;
		}
		return -1;
	}
}

int PairedEndReadsMerger::checkAlignment(string subjectString, string queryString, UINT16 rightminimumoverlap, UINT16 overallminimumoverlap,double tolerateErrorRate, bool direction)
{

	if(direction == true)
	{
		//check alignment with query 2, with query 1 hit by hash key
		//query 1: XXXXXXXXXXXKKKKKK query 2: MMMMMMMMXXXXXXXXX
		//subject:            KKKKKKXXXXXXXXXXMMMMMMMM
		if(queryString.length()<rightminimumoverlap||subjectString.length()<rightminimumoverlap) return -1;
		for(int i=subjectString.length()-rightminimumoverlap;i>=0;i--)
		{
			int overlap = subjectString.length()-i;
			int tolerateErrors = overlap*tolerateErrorRate;
			int errorCount = 0;
			for(int j=i;j<(int)subjectString.length() && j-i<(int)queryString.length()&& errorCount<=tolerateErrors;j++)
			{

				if(subjectString.at(j)!=queryString.at(j-i))errorCount++;

			}

			if(errorCount<=tolerateErrors) return i;
		}
		return -1;
	}
	else
	{
		//check alignment with query 2 (reverse), with query 1 (reverse) hit by hash key
		//query 2 (rev): XXXXXXXXXXXMMMMMM query 1 (rev): KKKKKKKKXXXXXXXXX
		//subject (rev):            MMMMMMXXXXXXXXXXXXXXXXKKKKKKKK
		if(queryString.length()<rightminimumoverlap||subjectString.length()<rightminimumoverlap) return -1;
		for(int i=rightminimumoverlap-1;i<(int)subjectString.length();i++)
		{
			int overlap = i+1;
			int tolerateErrors = overlap*tolerateErrorRate;
			int errorCount = 0;
			for(int j=i;j>=0&&i-j<(int)queryString.length()&& errorCount<=tolerateErrors;j--)
			{

				if(subjectString.at(j)!=queryString.at(queryString.length()-i-1+j))errorCount++;

			}

			if(errorCount<=tolerateErrors) return i;
		}
		return -1;
	}
}

bool PairedEndReadsMerger::start()
{
	if(Config::useErrorRate==true)
		return this->start_errorrate_seqfree();
	else
		return this->start_errornum_seqfree();
	return false;
}

//different implementation for start_errornum to save memory
bool PairedEndReadsMerger::start_errornum_seqfree()
{

	CLOCKSTART;
	MEMORYSTART;
	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * sRead = subjectReadList->at(i);
			SubjectEdge * sEdge = new SubjectEdge(sRead);
			subjectEdgeList->push_back(sEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTable(subjectEdge, true); //need a full scan here

			//clear the unused subject reads
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
//				cout<<subjectEdge->subjectRead->getName()<<"  "<<subjectEdge->alignmentList->size()<<endl;
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 int pairedEndAlign = -1;
					 QueryRead * rightQueryRead = alignment->queryRead->getPairEndRead()->rightRead;

					 if(alignment->queryOrientation==true)
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->getSequence(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength, this->maxMismatch, true);
					 }
					 else
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->reverseComplement(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength,this->maxMismatch,false);
					 }

					 if(pairedEndAlign!=-1)
					 {
						 int leftqueryend,subjectstart,subjectend, rightquerystart, rightqueryend;
						 bool suborient;
						 if(alignment->queryOrientation==true)
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = alignment->subjectStart;
							 subjectend = alignment->subjectEnd;
							 rightquerystart = pairedEndAlign+subjectstart;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = true;
						 }
						 else
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = -(alignment->subjectEnd)+alignment->queryEnd;
							 subjectend = -(alignment->subjectStart)+alignment->queryEnd;
							 rightquerystart = -(pairedEndAlign+alignment->subjectStart) + alignment->queryEnd;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = false;
						 }
						 if(leftqueryend-subjectstart+1+subjectend-rightquerystart+1>=this->overall_minimumOverlapLength)
						 {
							 PairedEndQueryAlignment* pairalign = new PairedEndQueryAlignment(alignment->queryRead->getPairEndRead(),alignment->subjectRead);
							 pairalign->left_queryEnd=leftqueryend;
							 pairalign->subjectStart=subjectstart;
							 pairalign->subjectEnd=subjectend;
							 pairalign->right_queryStart=rightquerystart;
							 pairalign->right_queryEnd=rightqueryend;
							 pairalign->subjectReadOrientation = suborient;
							#pragma omp critical(addAlignmentToQueryPair)
							 {
							 alignment->queryRead->getPairEndRead()->addToPairAlignmentList(pairalign);
							 }
						 }
					 }


					 delete alignment;//this alignment need to be deleted.

				}
			#pragma omp critical(AddSubjectRead)
				{
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
				}

			}
			delete subjectEdge;
			subjectEdgeList->at(i)=NULL;
		}

	}
		// end of omp parallel


	//need to process the internal sequence immediately per end reads pair
	vector<PairedEndRead*>* pairedReadList = this->hashTableMethod->getDataset()->getPairedEndReadList();
	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
	#pragma omp parallel
		{
	#pragma omp for schedule(dynamic)
			for(int i=0; i<(int)pairedReadList->size();i++)
			{
				PairedEndRead* pairedEndRead = pairedReadList->at(i);
				pairedEndRead->insertgapstringfrequencymap();
				pairedEndRead->clearPairEndAlignmentList();
			}
		}


		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

		subjectDataset->clearSubjectReadList();//it's done with all the subject reads in this chunk.

	}// end of while chunk loop


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
			pairedEndRead->reconstructSequencefromGapSeq(MINIMUMDEPTH, RECALLPERCENT);
		}
	}

	string notMergedFileName = Config::outputfilename+"notmerged.fasta";
	string mergedFileName = Config::outputfilename+"merged.fasta";
	printMergedToFile(notMergedFileName, mergedFileName);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;

	return true;

}

bool PairedEndReadsMerger::start_errorrate_seqfree()
{

	CLOCKSTART;
	MEMORYSTART;
	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * sRead = subjectReadList->at(i);
			SubjectEdge * sEdge = new SubjectEdge(sRead);
			subjectEdgeList->push_back(sEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTable(subjectEdge, true); //need a full scan here

			//clear the unused subject reads
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 int pairedEndAlign = -1;
					 QueryRead * rightQueryRead = alignment->queryRead->getPairEndRead()->rightRead;

					 if(alignment->queryOrientation==true)
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->getSequence(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength, this->maxErrorRate, true);
					 }
					 else
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->reverseComplement(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength,this->maxErrorRate,false);
					 }
					 if(pairedEndAlign!=-1)
					 {
						 int leftqueryend,subjectstart,subjectend, rightquerystart, rightqueryend;
						 bool suborient;
						 if(alignment->queryOrientation==true)
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = alignment->subjectStart;
							 subjectend = alignment->subjectEnd;
							 rightquerystart = pairedEndAlign+subjectstart;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = true;
						 }
						 else
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = -(alignment->subjectEnd)+alignment->queryEnd;
							 subjectend = -(alignment->subjectStart)+alignment->queryEnd;
							 rightquerystart = -(pairedEndAlign+alignment->subjectStart) + alignment->queryEnd;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = false;
						 }
						 if(leftqueryend-subjectstart+1+subjectend-rightquerystart+1>=this->overall_minimumOverlapLength)
						 {
							 PairedEndQueryAlignment* pairalign = new PairedEndQueryAlignment(alignment->queryRead->getPairEndRead(),alignment->subjectRead);
							 pairalign->left_queryEnd=leftqueryend;
							 pairalign->subjectStart=subjectstart;
							 pairalign->subjectEnd=subjectend;
							 pairalign->right_queryStart=rightquerystart;
							 pairalign->right_queryEnd=rightqueryend;
							 pairalign->subjectReadOrientation = suborient;
							#pragma omp critical(addAlignmentToQueryPair)
							 {
							 alignment->queryRead->getPairEndRead()->addToPairAlignmentList(pairalign);
							 }
						 }
					 }


					 delete alignment;//this alignment need to be deleted.

				}
			#pragma omp critical(AddSubjectRead)
				{
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
				}

			}
			delete subjectEdge;
			subjectEdgeList->at(i)=NULL;
		}

	}
		// end of omp parallel


	//need to process the internal sequence immediately per end reads pair
	vector<PairedEndRead*>* pairedReadList = this->hashTableMethod->getDataset()->getPairedEndReadList();
	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
	#pragma omp parallel
		{
	#pragma omp for schedule(dynamic)
			for(int i=0; i<(int)pairedReadList->size();i++)
			{
				PairedEndRead* pairedEndRead = pairedReadList->at(i);
				pairedEndRead->insertgapstringfrequencymap();
				pairedEndRead->clearPairEndAlignmentList();
			}
		}


		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

		subjectDataset->clearSubjectReadList();//it's done with all the subject reads in this chunk.

	}// end of while chunk loop


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
			pairedEndRead->reconstructSequencefromGapSeq(MINIMUMDEPTH, RECALLPERCENT);
		}
	}

	string notMergedFileName = Config::outputfilename+"notmerged.fasta";
	string mergedFileName = Config::outputfilename+"merged.fasta";
	printMergedToFile(notMergedFileName, mergedFileName);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;

	return true;

}

bool PairedEndReadsMerger::start_errornum()
{

	CLOCKSTART;
	MEMORYSTART;
	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * sRead = subjectReadList->at(i);
			SubjectEdge * sEdge = new SubjectEdge(sRead);
			subjectEdgeList->push_back(sEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTable(subjectEdge, true); //need a full scan here

			//clear the unused subject reads
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
//				cout<<subjectEdge->subjectRead->getName()<<"  "<<subjectEdge->alignmentList->size()<<endl;
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 int pairedEndAlign = -1;
					 QueryRead * rightQueryRead = alignment->queryRead->getPairEndRead()->rightRead;

					 if(alignment->queryOrientation==true)
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->getSequence(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength, this->maxMismatch, true);
					 }
					 else
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->reverseComplement(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength,this->maxMismatch,false);
					 }

					 if(pairedEndAlign!=-1)
					 {
						 int leftqueryend,subjectstart,subjectend, rightquerystart, rightqueryend;
						 bool suborient;
						 if(alignment->queryOrientation==true)
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = alignment->subjectStart;
							 subjectend = alignment->subjectEnd;
							 rightquerystart = pairedEndAlign+subjectstart;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = true;
						 }
						 else
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = -(alignment->subjectEnd)+alignment->queryEnd;
							 subjectend = -(alignment->subjectStart)+alignment->queryEnd;
							 rightquerystart = -(pairedEndAlign+alignment->subjectStart) + alignment->queryEnd;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = false;
						 }
						 if(leftqueryend-subjectstart+1+subjectend-rightquerystart+1>=this->overall_minimumOverlapLength)
						 {
							 PairedEndQueryAlignment* pairalign = new PairedEndQueryAlignment(alignment->queryRead->getPairEndRead(),alignment->subjectRead);
							 pairalign->left_queryEnd=leftqueryend;
							 pairalign->subjectStart=subjectstart;
							 pairalign->subjectEnd=subjectend;
							 pairalign->right_queryStart=rightquerystart;
							 pairalign->right_queryEnd=rightqueryend;
							 pairalign->subjectReadOrientation = suborient;
							#pragma omp critical(addAlignmentToQueryPair)
							 {
							 alignment->queryRead->getPairEndRead()->addToPairAlignmentList(pairalign);
							 }
						 }
					 }


					 delete alignment;//this alignment need to be deleted.

				}
			#pragma omp critical(AddSubjectRead)
				{
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
				}

			}
			delete subjectEdge;
			subjectEdgeList->at(i)=NULL;
		}

	}
		// end of omp parallel


		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

	}// end of while chunk loop


	vector<PairedEndRead*>* pairedReadList = this->hashTableMethod->getDataset()->getPairedEndReadList();

omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(int i=0; i<(int)pairedReadList->size();i++)
		{
			PairedEndRead* pairedEndRead = pairedReadList->at(i);
			//int minimumDepth = 5;
			if(pairedEndRead->insertMatrixToMap(MINIMUMDEPTH))
				pairedEndRead->reconstructSequence(MINIMUMDEPTH, RECALLPERCENT);

		}
	}

	string notMergedFileName = Config::outputfilename+"notmerged.fasta";
	string mergedFileName = Config::outputfilename+"merged.fasta";
	printMergedToFile(notMergedFileName, mergedFileName);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;

	return true;

}


bool PairedEndReadsMerger::start_errorrate()
{

	CLOCKSTART;
	MEMORYSTART;
	vector<SubjectRead *>* subjectReadList= new vector<SubjectRead *>();
	subjectReadList->clear();

	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	while(subjectDataset->loadNextChunk(subjectReadList))
	{
		vector<SubjectEdge*>* subjectEdgeList = new  vector<SubjectEdge*>();
		for(UINT64 i = 0; i < subjectReadList->size(); i++)
		{
			SubjectRead * sRead = subjectReadList->at(i);
			SubjectEdge * sEdge = new SubjectEdge(sRead);
			subjectEdgeList->push_back(sEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTable(subjectEdge, true); //need a full scan here

			//clear the unused subject reads
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
//				cout<<subjectEdge->subjectRead->getName()<<"  "<<subjectEdge->alignmentList->size()<<endl;
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 int pairedEndAlign = -1;
					 QueryRead * rightQueryRead = alignment->queryRead->getPairEndRead()->rightRead;

					 if(alignment->queryOrientation==true)
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->getSequence(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength, this->maxErrorRate, true);
					 }
					 else
					 {
						 pairedEndAlign = this->checkAlignment(subjectEdge->subjectRead->getSequence(), rightQueryRead->reverseComplement(),this->right_minimumOverlapLength, this->overall_minimumOverlapLength,this->maxErrorRate,false);
					 }

					 if(pairedEndAlign!=-1)
					 {
						 int leftqueryend,subjectstart,subjectend, rightquerystart, rightqueryend;
						 bool suborient;
						 if(alignment->queryOrientation==true)
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = alignment->subjectStart;
							 subjectend = alignment->subjectEnd;
							 rightquerystart = pairedEndAlign+subjectstart;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = true;
						 }
						 else
						 {
							 leftqueryend = alignment->queryEnd;
							 subjectstart = -(alignment->subjectEnd)+alignment->queryEnd;
							 subjectend = -(alignment->subjectStart)+alignment->queryEnd;
							 rightquerystart = -(pairedEndAlign+alignment->subjectStart) + alignment->queryEnd;
							 rightqueryend = rightquerystart-1+rightQueryRead->getReadLength();
							 suborient = false;
						 }
						 if(leftqueryend-subjectstart+1+subjectend-rightquerystart+1>=this->overall_minimumOverlapLength)
						 {
							 PairedEndQueryAlignment* pairalign = new PairedEndQueryAlignment(alignment->queryRead->getPairEndRead(),alignment->subjectRead);
							 pairalign->left_queryEnd=leftqueryend;
							 pairalign->subjectStart=subjectstart;
							 pairalign->subjectEnd=subjectend;
							 pairalign->right_queryStart=rightquerystart;
							 pairalign->right_queryEnd=rightqueryend;
							 pairalign->subjectReadOrientation = suborient;
							#pragma omp critical(addAlignmentToQueryPair)
							 {
							 alignment->queryRead->getPairEndRead()->addToPairAlignmentList(pairalign);
							 }
						 }
					 }


					 delete alignment;//this alignment need to be deleted.

				}
			#pragma omp critical(AddSubjectRead)
				{
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
				}

			}
			delete subjectEdge;
			subjectEdgeList->at(i)=NULL;
		}

	}
		// end of omp parallel


		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

	}// end of while chunk loop


	vector<PairedEndRead*>* pairedReadList = this->hashTableMethod->getDataset()->getPairedEndReadList();

omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(int i=0; i<(int)pairedReadList->size();i++)
		{
			PairedEndRead* pairedEndRead = pairedReadList->at(i);
			//int minimumDepth = 5;
			if(pairedEndRead->insertMatrixToMap(MINIMUMDEPTH))
				pairedEndRead->reconstructSequence(MINIMUMDEPTH, RECALLPERCENT);
		}
	}

	string notMergedFileName = Config::outputfilename+"notmerged.fasta";
	string mergedFileName = Config::outputfilename+"merged.fasta";
	printMergedToFile(notMergedFileName, mergedFileName);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;

	return true;

}

bool PairedEndReadsMerger::printMergedToFile(string notMergedFileName, string mergedFileName)
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



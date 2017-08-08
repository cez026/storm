/*
 * OverlapGraphConstructor.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#include "OverlapGraphConstructor.h"

OverlapGraphConstructor::OverlapGraphConstructor(HashTableMethod * hashTableMethod) {
	// TODO Auto-generated constructor stub
	this->totaledgenumber = 0;
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashTableMethod = hashTableMethod;

}

OverlapGraphConstructor::~OverlapGraphConstructor() {
	// TODO Auto-generated destructor stub

}


bool OverlapGraphConstructor::start()
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
			this->hashTableMethod->searchHashTable(subjectEdge);

			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				bool newalignadded=false;
				//add alignments to the query reads before the edges are destroyed.
				for(INT64 j = 0; j < (int)subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 if(!this->isContainedAlignment(alignment))//only add non-contained read;
					 {
						 newalignadded = true;
					#pragma omp critical(addAlignmentToQueryRead)
					 {
					 alignment->queryRead->addAlignment(alignment);
					 this->totaledgenumber++;
					 }
					 }
					 else
					 {
						 delete alignment;
					 }

				}
				if(newalignadded == true)
				{
				#pragma omp critical(AddSubjectRead)
				{
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
				}
				}
//				 cout<<"add subject:"<<subjectEdge->subjectRead->getName()<<endl;
			}
			delete subjectEdge;
			subjectEdgeList->at(i)=NULL;
		}

	}
		// end of omp parallel



		//
		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

	}// end of while chunk loop

	printEdgesToFile(Config::outputfilename);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;
cout<<"edge: "<<this->totaledgenumber<<endl;

	return true;

}

bool OverlapGraphConstructor::start(bool transitiveEdgeRemoval)
{

	CLOCKSTART;
	MEMORYSTART;


	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	subjectDataset->setReadFilterSize(Config::minimumOverlapLength);
	vector<SubjectEdge*>* subjectEdgesList = new  vector<SubjectEdge*>();
	subjectEdgesList->clear();
	while(subjectDataset->loadNextChunk(subjectEdgesList))
	{
//		{
//			CURMEM;
//			cout<<"searching is started"<<endl;
//			CURTIME;
//		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgesList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgesList->at(i);
			//use full pair-wise comparison, redundancy is removed by only considering one side alignment in this case.
			this->hashTableMethod->searchHashTableFullScan(subjectEdge);

			if(subjectEdge->alignmentList!=NULL)
			{
				INT64 usefullinks=0;
				for(INT64 j = 0; j < (int)subjectEdge->alignmentList->size(); j++)
				{
					Alignment * alignment= subjectEdge->alignmentList->at(j);
					QueryRead * queryreadpointer = alignment->queryRead;
					if(!this->isContainedAlignment(alignment))//only add non-contained read;
					{
						queryreadpointer->addAlignmentWithLock(alignment);
//#pragma omp critical(addalignment)
//						{
//							queryreadpointer->addAlignment(alignment);
//						}
						usefullinks++;
//						alignment->subjectRead->increaselinknumberWithLock();
					}
					else
					{
						delete alignment;
					}
					 //after using this alignment, set it to null
					 subjectEdge->alignmentList->at(j)=NULL;
				}
				subjectEdge->subjectRead->setLinkNum(usefullinks);
			}
			//clear the alignment list right after processing, to save memory
			subjectEdge->clearAlignmentList();

		}//end of for loop

	}// end of omp parallel


//	{
//		CURMEM;
//		cout<<"searching is done"<<endl;
//		CURTIME;
//	}
	//--------transitive edge removal-----------
	if(transitiveEdgeRemoval==true)
	{
	//need to process transitive edge removal immediately per query read
	QueryDataset * dataset = this->hashTableMethod->getDataset();
	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
	#pragma omp parallel
		{
	#pragma omp for schedule(dynamic)
			for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
			{

				QueryRead* queryRead = dataset->getReadFromID(i);
//#pragma omp critical(printedge)
//				{
//				if(i<100&&queryRead->getAlignmentNumber()>0)cout<<queryRead->getIdentifier()<<":"<<queryRead->getAlignmentNumber()<<endl;
//				}
				queryRead->markTransitiveEdgeIncremental();
				queryRead->deleteTransitiveRemovableRedundantEdge();
//#pragma omp critical(printedge)
//				{
//				if(i<100&&queryRead->getAlignmentNumber()>0)cout<<queryRead->getIdentifier()<<":"<<queryRead->getAlignmentNumber()<<endl;
//				}
			}
		}//end of parallel
//		{
//			CURMEM;
//			cout<<"transitive removal is done"<<endl;
//			CURTIME;
//		}
		//need to delete un-needed subject read in subject dataset
		subjectDataset->clearNoLinkageReadsAndDelete();
//		{
//			CURMEM;
//			cout<<"subject database removal is done"<<endl;
//			CURTIME;
//		}
	}//end of if



	//--------add necessary subject read, and delete rest of those-------
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
	for(UINT64 i = 0; i < subjectEdgesList->size(); i++)
	{
		SubjectEdge * subjectEdge = subjectEdgesList->at(i);
		if(subjectEdge->subjectRead->isNoLinkage())
		{
			delete subjectEdge->subjectRead;
		}
		else
		{
#pragma omp critical(AddSubjectRead)
			{
			subjectDataset->addSubjectRead(subjectEdge->subjectRead);
			}
		}
		subjectEdge->subjectRead = NULL;
	}
	}//end of parallel
	//clear the subject read vector for the next round
	subjectEdgesList->clear();
	subjectEdgesList->resize(0);
//	{
//		CURMEM;
//		cout<<"clear is done"<<endl;
//		CURTIME;
//		cout<<"*****one round is done*******"<<endl;
//	}
	}// end of while chunk loop


	delete subjectEdgesList;

	printEdgesToFileWithoutDuplicate(Config::outputfilename);

	delete subjectDataset;

	MEMORYSTOP;
	CLOCKSTOP;

	return true;

}



bool OverlapGraphConstructor::printEdgesToFile(string outFileName)
{
	CLOCKSTART;
	MEMORYSTART;
	ofstream filePointer;
	filePointer.open(outFileName.c_str());
	if(filePointer == NULL)
	{
		cout<<"Unable to open file: "<<endl;
		return false;
	}
	QueryDataset * dataset = this->hashTableMethod->getDataset();


	if(Config::numberOfThreads==1)
	{
		for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
			dataset->getReadFromID(i)->printAlignmentToFile(filePointer);
	}
	else
	{
		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel shared(filePointer)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
				{
					string outputstring = dataset->getReadFromID(i)->getStringForPrintAlignmentToFile();
				#pragma omp critical(printedge)
					{
					filePointer<<outputstring;
					}
				}
			}
	}



	filePointer.close();
	MEMORYSTOP;
	CLOCKSTOP;
	return true;
}

bool OverlapGraphConstructor::printEdgesToFileWithoutDuplicate(string outFileName)
{
	CLOCKSTART;
	MEMORYSTART;
	ofstream filePointer;
	filePointer.open(outFileName.c_str());
	if(filePointer == NULL)
	{
		cout<<"Unable to open file: "<<endl;
		return false;
	}
	QueryDataset * dataset = this->hashTableMethod->getDataset();



		omp_set_dynamic(0);
		omp_set_num_threads(Config::numberOfThreads);
		#pragma omp parallel shared(filePointer)
			{
		#pragma omp for schedule(dynamic)
				for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
				{
					string outputstring = dataset->getReadFromID(i)->getStringForPrintAlignmentToFileWithoutDuplicate();
				#pragma omp critical(printedge)
					{
					filePointer<<outputstring;
					}
				}
			}


	filePointer.close();
	MEMORYSTOP;
	CLOCKSTOP;
	return true;
}


bool OverlapGraphConstructor::printEdgesToFileWithoutTransitiveEdge(string outFileName)
{
	CLOCKSTART;
	MEMORYSTART;
	ofstream filePointer;
	filePointer.open(outFileName.c_str());
	if(filePointer == NULL)
	{
		cout<<"Unable to open file: "<<endl;
		return false;
	}
	QueryDataset * dataset = this->hashTableMethod->getDataset();

	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
	#pragma omp parallel shared(filePointer)
		{
	#pragma omp for schedule(dynamic)
			for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
			{
				string outputstring = dataset->getReadFromID(i)->getStringForPrintNonTransitiveAlignmentToFile();
			#pragma omp critical(printedge)
				{
//					cout<<outputstring<<" "<<i<<endl;
				filePointer<<outputstring;
				}
			}
		}


	filePointer.close();
	MEMORYSTOP;
	CLOCKSTOP;
	return true;
}

//either subject is contained in query or query is contained in subject
bool OverlapGraphConstructor::isContainedAlignment(Alignment * subjectAlignment)
{
	if(subjectAlignment->subjectStart<=0&&subjectAlignment->subjectEnd>=subjectAlignment->queryEnd)
		return true;
	else if(subjectAlignment->subjectStart>=0&&subjectAlignment->subjectEnd<=subjectAlignment->queryEnd)
		return true;
	else return false;
}

/*
 *

bool OverlapGraphConstructor::start(bool transitiveEdgeRemoval)
{

	CLOCKSTART;
	MEMORYSTART;


	SubjectDataset *subjectDataset = new SubjectDataset();
	subjectDataset->setFilenameList(Config::subjectFilenameList);
	subjectDataset->setReadFilterSize(Config::minimumOverlapLength);
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
			//use full pair-wise comparison, redundancy is removed by only considering one side alignment in this case.
			this->hashTableMethod->searchHashTableFullScan(subjectEdge);

			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				bool newalignadded=false;
				//add alignments to the query reads before the edges are destroyed.
//#pragma omp critical(print)
//				{
//				cout<<subjectEdge->alignmentList->size()<<"  "<<subjectEdge->subjectRead->getName()<<endl;
//				}

				for(INT64 j = 0; j < (int)subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);


					 if(!this->isContainedAlignment(alignment))//only add non-contained read;
					 {

					#pragma omp critical(addAlignmentToQueryRead)
					 {
						 alignment->queryRead->addAlignment(alignment);
						 this->totaledgenumber++;
						 if(newalignadded==false)
						 {
							 subjectDataset->addSubjectRead(alignment->subjectRead);
							 newalignadded = true;
						 }
					 }
					 alignment->subjectRead->increaselinknumber();

					 }
					 else
					 {
						 delete alignment;
					 }
					 subjectEdge->alignmentList->at(j)=NULL;

				}
//				if(newalignadded == true)
//				{
//				#pragma omp critical(AddSubjectRead)
//				{
//				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
//				}
//				}


//				 cout<<"add subject:"<<subjectEdge->subjectRead->getName()<<endl;

			}
			delete subjectEdge;
			subjectEdgesList->at(i)=NULL;
		}

	}
		// end of omp parallel

	//clear the subject read vector for the next round
	subjectEdgesList->clear();
	subjectEdgesList->resize(0);

	if(transitiveEdgeRemoval==true)
	{
	//need to process transitive edge removal immediately per query read
	QueryDataset * dataset = this->hashTableMethod->getDataset();
	omp_set_dynamic(0);
	omp_set_num_threads(Config::numberOfThreads);
	#pragma omp parallel
		{
	#pragma omp for schedule(dynamic)
			for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
			{
				QueryRead* queryRead = dataset->getReadFromID(i);
				queryRead->markTransitiveEdge();
				queryRead->deleteTransitiveRemovableRedundantEdge();
			}
		}

		subjectDataset->clearNoLinkageReads();//try to clear subject reads as much as possible to save memory

	}//end of if
cout<<"end of if"<<endl;
	}// end of while chunk loop

	delete subjectDataset;
	delete subjectEdgesList;

	printEdgesToFileWithoutDuplicate(Config::outputfilename);


	MEMORYSTOP;
	CLOCKSTOP;
cout<<"edge: "<<this->totaledgenumber<<endl;

	return true;

}

//---------------------------------------------------------------------------------

bool OverlapGraphConstructor::start(bool transitiveEdgeRemoval)
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
			//use full pair-wise comparison, redundancy is removed by only considering one side alignment in this case.
			this->hashTableMethod->searchHashTableFullScan(subjectEdge);

			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				bool newalignadded=false;
				//add alignments to the query reads before the edges are destroyed.
				for(INT64 j = 0; j < subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 if(!this->isContainedAlignment(alignment))//only add non-contained read;
					 {
						 newalignadded = true;
					#pragma omp critical(addAlignmentToQueryRead)
					 {
					 alignment->queryRead->addAlignment(alignment);
					 this->totaledgenumber++;
					 }
					 }
					 else
					 {
						 delete alignment;
					 }

				}
				if(newalignadded == true)
				{
				#pragma omp critical(AddSubjectRead)
				{
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
				}
				}
//				 cout<<"add subject:"<<subjectEdge->subjectRead->getName()<<endl;
			}
			delete subjectEdge;
			subjectEdgeList->at(i)=NULL;
		}

	}
		// end of omp parallel



		//
		subjectEdgeList->clear();
		delete subjectEdgeList;
		//clear the subject read vector for the next round
		subjectReadList->clear();
		subjectReadList->resize(0);

	}// end of while chunk loop

	if(transitiveEdgeRemoval==false)
	printEdgesToFileWithoutDuplicate(Config::outputfilename);
	else
	printEdgesToFileWithoutTransitiveEdge(Config::outputfilename);
	delete subjectReadList;
	delete subjectDataset;
	MEMORYSTOP;
	CLOCKSTOP;
cout<<"edge: "<<this->totaledgenumber<<endl;

	return true;

}
 */

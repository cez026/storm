/*
 * GeneralQueryDatasetFilter.cpp
 *
 *  Created on: Apr 29, 2015
 *      Author: qy2
 */

#include "GeneralQueryDatasetFilter.h"

GeneralQueryDatasetFilter::GeneralQueryDatasetFilter(HashTableMethod * hashTableMethod) {
	// TODO Auto-generated constructor stub

	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashTableMethod = hashTableMethod;
}

GeneralQueryDatasetFilter::~GeneralQueryDatasetFilter() {
	// TODO Auto-generated destructor stub
}

bool GeneralQueryDatasetFilter::start()
{
	string outFileName = Config::outputfilename;
	filePointer.open(outFileName.c_str());
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
			this->hashTableMethod->searchHashTableFullScanAndMarkQueryRead(subjectEdge); //and also mark the contained read in the query read structure
			delete subjectEdge->subjectRead;
			delete subjectEdge;
			subjectEdgesList->at(i) = NULL;
		}

	}
		// end of omp parallel

		//

		//clear the subject read vector for the next round
		subjectEdgesList->clear();
		subjectEdgesList->resize(0);

	}// end of while chunk loop
	this->printToFile();

	delete subjectEdgesList;
	delete subjectDataset;

	MEMORYSTOP;
	CLOCKSTOP;

filePointer.close();
	return true;
}

void GeneralQueryDatasetFilter::printToFile()
{
	CLOCKSTART;
	MEMORYSTART;
	UINT64 count = 0;
	QueryDataset * dataset = this->hashTableMethod->getDataset();
	if(Config::useID == false)
	{
		for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
		{
			QueryRead *queryRead = dataset->getReadFromID(i);
			if(queryRead->flag4Removal==false)
			{
				count++;
				std::stringstream sstm;
				sstm<<">"<<queryRead->getName()<<endl;
				sstm<<queryRead->getSequence()<<endl;
				this->filePointer<<sstm.str();
				if(count%100000==0)	// Show the progress.
					cout<<"counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;
			}
		}
	}
	else
	{
		if(Config::isPrintIDmap==false)
		{
			for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
			{
				QueryRead *queryRead = dataset->getReadFromID(i);
				if(queryRead->flag4Removal==false)
				{
					count++;
					std::stringstream sstm;
					UINT64 queryReadID = queryRead->getIdentifier()+Config::startID;
					sstm<<">"<<queryReadID<<endl;
					sstm<<queryRead->getSequence()<<endl;
					this->filePointer<<sstm.str();
					if(count%100000==0)	// Show the progress.
						cout<<"counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;
				}
			}
		}
		else
		{
			string mapFileName = Config::IDmapfilename;
			ofstream mapfilePointer;
			mapfilePointer.open(mapFileName.c_str());
			for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
			{
				QueryRead *queryRead = dataset->getReadFromID(i);
				if(queryRead->flag4Removal==false)
				{
					count++;
					std::stringstream sstm;
					UINT64 queryReadID = queryRead->getIdentifier()+Config::startID;
					sstm<<">"<<queryReadID<<endl;
					sstm<<queryRead->getSequence()<<endl;
					this->filePointer<<sstm.str();
					if(count%100000==0)	// Show the progress.
						cout<<"counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;

					std::stringstream sstm2;
					sstm2<<queryReadID<<"\t"<<queryRead->getName()<<endl;
					mapfilePointer<<sstm2.str();
				}

			}
			cout<<"ID mapping file is generated"<<endl;
			mapfilePointer.close();
		}
	}
	cout<<"Total counter: " << setw(10) << count << " Reads Written. " << setw(10) << endl;
	MEMORYSTOP;
	CLOCKSTOP;
}
/*
bool GeneralQueryDatasetFilter::start()
{
	string outFileName = Config::outputfilename;
	filePointer.open(outFileName.c_str());
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
			SubjectRead * subjectRead = subjectReadList->at(i);
			SubjectEdge * subjectEdge = new SubjectEdge(subjectRead);
			subjectEdgeList->push_back(subjectEdge);
		}
omp_set_dynamic(0);
omp_set_num_threads(Config::numberOfThreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			this->hashTableMethod->searchHashTableFullScanAndMarkQueryRead(subjectEdge); //and also mark the contained read in the query read structure
			delete subjectEdge->subjectRead;
			delete subjectEdge;
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
	this->printToFile();

	delete subjectReadList;
	delete subjectDataset;

	MEMORYSTOP;
	CLOCKSTOP;

filePointer.close();
	return true;
}
*/

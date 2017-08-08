/*
 * OmegaGraphConstructor.cpp
 *
 *  Created on: Feb 16, 2015
 *      Author: qy2
 */

#include "OmegaGraphConstructor.h"

OmegaGraphConstructor::OmegaGraphConstructor(OmegaHashTable * omegaHashTable)
{
	// TODO Auto-generated constructor stub
	this->omegaHashTable = omegaHashTable;
	this->totaledgenumber = 0;
	this->minimumOverlapLength = Config::minimumOverlapLength;
}

OmegaGraphConstructor::~OmegaGraphConstructor() {
	// TODO Auto-generated destructor stub
}
/*
bool OmegaGraphConstructor::start()
{
	this->totaledgenumber = 0;
	filePointer.open("omegatemp.txt");
	CLOCKSTART;
	MEMORYSTART;
//	estimatedGenomeSize = 0;
//	numberOfNodes = 0;
//	numberOfEdges = 0;
//	flowComputed = false;
	OmegaHashTable * hashTable = this->omegaHashTable;
	QueryDataset * dataSet = hashTable->getDataset();
	UINT64 counter = 0;
	vector<int> *exploredReads = new vector<int>; // 0--unexplored
	exploredReads->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<UINT64> * queue = new vector<UINT64>;
	queue->reserve(dataSet->getNumberOfUniqueReads()+1);

//	vector<markType> *markedNodes = new vector<markType>;
//	markedNodes->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector< vector<Edge *> * >* graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);



	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // Initialization
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
		exploredReads->push_back(0);
		queue->push_back(0);
//		markedNodes->push_back(VACANT);
	}


	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(exploredReads->at(i) == 0)
		{
			UINT64 start = 0, end = 0; 											// Initialize queue start and end.
			queue->at(end++) = i;
			while(start < end) 													// This loop will explore all connected component starting from read i.
			{
				counter++;
				UINT64 read1 = queue->at(start++);
				if(exploredReads->at(read1) == 0)
				{
					insertAllEdgesOfRead(read1, exploredReads);					// Explore current node.
					exploredReads->at(read1) = 1;
				}
				if(graph->at(read1)->size() != 0) 								// Read has some edges (required only for the first read when a new queue starts.
				{
					if(exploredReads->at(read1) == 1) 					// Explore unexplored neighbors first.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++ )
						{
							UINT64 read2 = graph->at(read1)->at(index1)->destination->getIdentifier();
							if(exploredReads->at(read2) == 0) 			// Not explored.
							{
								queue->at(end++) = read2; 						// Put in the queue.
								insertAllEdgesOfRead(read2, exploredReads);
								exploredReads->at(read2) = 1;
							}
						}
//						exploredReads->at(read1) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
					}
				}
				if(counter%100000==0)	// Show the progress.
					cout<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << endl;
			}
		}
	}
	cout<<"total counter: " << setw(10) << counter << " Nodes: " << setw(10)  << endl;
	cout<<"total edge number: "<<this->totaledgenumber<<endl;
	delete exploredReads;
	delete queue;
//	delete markedNodes;
	filePointer.close();

	MEMORYSTOP;
	CLOCKSTOP;
	return true;
}


bool OmegaGraphConstructor::insertAllEdgesOfRead(UINT64 readNumber, vector<int> * exploredReads)
{
	QueryRead *read1 = this->omegaHashTable->getDataset()->getReadFromID(readNumber); 	// Get the current read read1.
	string readString = read1->getSequence(); 		// Get the forward string of read1.
	string subString;
//YAO	for(UINT64 j = 1; j < read1->getReadLength()-hashTable->getHashStringLength(); j++) // For each proper substring of length getHashStringLength of read1

	for(UINT64 j = 0; j <= read1->getReadLength()-this->omegaHashTable->getHashStringLength(); j++) // For each proper substring of length getHashStringLength of read1
	{
		subString = readString.substr(j,this->omegaHashTable->getHashStringLength());  // Get the proper substring s of read1.
		vector<UINT64> * listOfReads=this->omegaHashTable->getListOfReads(subString); // Search the string in the hash table.

		if(!listOfReads->empty()) // If there are some reads that contain s as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT16 overlapOffset;
				UINT8 orientation;
				QueryRead *read2 = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(exploredReads->at(read2->getIdentifier())!= 0) 			// No need to discover the same edge again. All edges of read2 is already inserted in the graph.
						continue;
				if(read1->flag4Removal==false && read2->flag4Removal==false && checkOverlap(read1,read2,(data >> 62),j)) // Both read need to be non contained.
				{
				std::stringstream sstm;

				int orientation = data>>62;
//				cout<<read1->getReadNumber() <<" "<< " >>> " << read2->getReadNumber() <<" orientation: "<<orientation<<" position: "<<j<<endl;

//				if(read1->getStringForward()<read2->getStringForward())
//				sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<read1->getStringForward()<<"||||"<<read2->getStringForward()<<endl;
//				else
//					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<read2->getStringForward()<<"||||"<<read1->getStringForward()<<endl;
//				filePointer<<sstm.str();
					this->totaledgenumber++;
				}
			}
		}
	}
//YAO	if(graph->at(readNumber)->size() != 0)
//YAO		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}
*/
/*
bool OmegaGraphConstructor::start()
{
	CLOCKSTART;
	MEMORYSTART;

	QueryDataset * querydataset = this->omegaHashTable->getDataset();
	for(UINT64 i=1;i<=querydataset->getNumberOfUniqueReads();i++)
	{
		QueryRead *curRead = querydataset->getReadFromID(i); 	// Get the current read subject read.
		string subjectReadString =curRead->getSequence(); 		// Get the forward string of subject read.
		string subString;
		for(UINT64 j = 0; j <= curRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
		{
			subString = subjectReadString.substr(j,this->omegaHashTable->getHashStringLength());
			vector<UINT64> * listOfReads=this->omegaHashTable->getListOfReads(subString); // Search the string in the hash table.
//			vector<UINT64> * listOfReads;
//			if(i<querydataset->getNumberOfUniqueReads()/2)
//				listOfReads=this->omegaHashTable->getListOfReads(subString);
//			else listOfReads=this->omegaHashTable->hashTable->at(0);
			if(!listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
			{
				for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
				{
					UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
					UINT16 overlapOffset;
					UINT8 orientation;
					QueryRead *qRead = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
					if(curRead->getName()<qRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						continue;
					}

					if(curRead->flag4Removal==false && qRead->flag4Removal==false && checkOverlap(curRead,qRead,(data >> 62),j)) // Both read need to be non contained.
					{
						this->totaledgenumber++;

					}
				}
			}
		}
	}

	MEMORYSTOP;
	CLOCKSTOP;
	cout<<"edge: "<<this->totaledgenumber<<endl;
	filePointer.close();
	return true;
}
*/

bool OmegaGraphConstructor::start()
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
			this->searchHashTable(subjectEdge);
		}

	}
		// end of omp parallel


	//clear the unused subject reads
		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				//add alignments to the query reads before the edges are destroyed.
				for(INT64 j = 0; j < (int)subjectEdge->alignmentList->size(); j++)
				{
					 Alignment * alignment= subjectEdge->alignmentList->at(j);
					 alignment->queryRead->addAlignment(alignment);
					 this->totaledgenumber++;

				}
				 subjectDataset->addSubjectRead(subjectEdge->subjectRead);
//				 cout<<"add subject:"<<subjectEdge->subjectRead->getName()<<endl;
			}
			delete subjectEdge;

		}

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

bool OmegaGraphConstructor::start(bool transitiveEdgeRemoval)
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
			this->searchHashTableFullScan(subjectEdge);

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
	QueryDataset * dataset = this->omegaHashTable->getDataset();
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




bool OmegaGraphConstructor::printEdgesToFile(string outFileName)
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
	QueryDataset * dataset = this->omegaHashTable->getDataset();
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

bool OmegaGraphConstructor::printEdgesToFileWithoutDuplicate(string outFileName)
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
	QueryDataset * dataset = this->omegaHashTable->getDataset();

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

bool OmegaGraphConstructor::printEdgesToFileWithoutTransitiveEdge(string outFileName)
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
	QueryDataset * dataset = this->omegaHashTable->getDataset();

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
				filePointer<<outputstring;
				}
			}
		}


	filePointer.close();
	MEMORYSTOP;
	CLOCKSTOP;
	return true;
}

bool OmegaGraphConstructor::searchHashTable(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
//	for(UINT64 j = 0; j <= subjectRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
	for(UINT64 j = 1; j < subjectRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
//avoid self non trivial alignment
	{
		subString = subjectReadString.substr(j,this->omegaHashTable->getHashStringLength());
		vector<UINT64> * listOfReads=this->omegaHashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
//				UINT16 overlapOffset;
//				UINT8 orientation;
				QueryRead *queryRead = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(subjectRead->getName()<queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
				{
//					cout << sRead->getName() << " <<<<<< " << qRead->getName() << " ????? " <<j<<" and " << k << endl;
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}

				if(subjectRead->flag4Removal==false && queryRead->flag4Removal==false && checkOverlap(subjectRead,queryRead,(data >> 62),j)) // Both read need to be non contained.
				{
				int orientation = data >> 62;
//YAO					cout<<sRead->getName() <<" "<< " >>> " << qRead->getName() <<" orientation: "<<orientation<<" position: "<<j<<endl;

//YAO					cout << "S  " << sRead->getSequence() << endl;
//YAO				cout << "Q  " << qRead->getSequence() << endl;
//YAO					std::stringstream sstm;
//YAO					if(sRead->getSequence()<qRead->getSequence())
//YAO					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<sRead->getSequence()<<"||||"<<qRead->getSequence()<<endl;
//YAO					else
//YAO					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<qRead->getSequence()<<"||||"<<sRead->getSequence()<<endl;
//YAO				filePointer<<sstm.str();

					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					switch (orientation) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
					case 0:
						subjectAlignment->queryOrientation = true;
						subjectAlignment->subjectStart = -j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					case 1:
						subjectAlignment->queryOrientation = true;
						subjectAlignment->subjectStart = queryRead->getReadLength()-this->omegaHashTable->getHashStringLength()-j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					case 2:
						subjectAlignment->queryOrientation = false;
						subjectAlignment->subjectStart = -j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					case 3:
						subjectAlignment->queryOrientation = false;
						subjectAlignment->subjectStart = queryRead->getReadLength()-this->omegaHashTable->getHashStringLength()-j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					default:;
					}
					subjectEdge->addAlignment(subjectAlignment);
				}
			}
		}
	}
//YAO	if(graph->at(readNumber)->size() != 0)
//YAO		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdges); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}

bool OmegaGraphConstructor::searchHashTableFullScan(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
//	for(UINT64 j = 0; j <= subjectRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
	for(UINT64 j = 1; j < subjectRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
//avoid self non trivial alignment
	{
		subString = subjectReadString.substr(j,this->omegaHashTable->getHashStringLength());
		vector<UINT64> * listOfReads=this->omegaHashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
//				UINT16 overlapOffset;
//				UINT8 orientation;
				QueryRead *queryRead = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(subjectRead->getName()==queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
				{
//					cout << sRead->getName() << " <<<<<< " << qRead->getName() << " ????? " <<j<<" and " << k << endl;
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}

				if(subjectRead->flag4Removal==false && queryRead->flag4Removal==false && checkOverlap(subjectRead,queryRead,(data >> 62),j)) // Both read need to be non contained.
				{
				int orientation = data >> 62;

					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					switch (orientation) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
					case 0:
						subjectAlignment->queryOrientation = true;
						subjectAlignment->subjectStart = -j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					case 1:
						subjectAlignment->queryOrientation = true;
						subjectAlignment->subjectStart = queryRead->getReadLength()-this->omegaHashTable->getHashStringLength()-j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					case 2:
						subjectAlignment->queryOrientation = false;
						subjectAlignment->subjectStart = -j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					case 3:
						subjectAlignment->queryOrientation = false;
						subjectAlignment->subjectStart = queryRead->getReadLength()-this->omegaHashTable->getHashStringLength()-j;
						subjectAlignment->queryEnd = queryRead->getReadLength()-1;
						subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
						break;
					default:;
					}
					subjectEdge->addAlignment(subjectAlignment);
				}
			}
		}
	}
	return true;
}

/**********************************************************************************************************************
	Checks if two read overlaps.
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	Orientation 0 means prefix of forward of the read2
	Orientation 1 means suffix of forward of the read2
	Orientation 2 means prefix of reverse of the read2
	Orientation 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read1 and read2 overlap.
**********************************************************************************************************************/
bool OmegaGraphConstructor::checkOverlap(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getSequence(); // Get the forward string of read1
	UINT64 hashStringLength = this->omegaHashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement(); // Get the string from read2 according to orient.
	if(orient == 0 || orient == 2)		// orient 0
										//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
										//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
										//				OR
										// orient 2
										//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
										//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(string1.length()-start<this->minimumOverlapLength)return false;
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
//		return string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)); // If the remaining strings match.
		return strcmp(string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)).c_str(), string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)).c_str())==0; // If the remaining strings match.

	}
	else								// orient 1
										//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
										//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
										//				OR
										// orient 3
										//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
										//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(start+hashStringLength<this->minimumOverlapLength)return false;
		if(string2.length()-hashStringLength < start)
			return false;
//		return string1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
		return strcmp(string1.substr(0, start).c_str(), string2.substr(string2.length()-hashStringLength-start, start).c_str())==0; // If the remaining strings match.

	}
}

bool OmegaGraphConstructor::checkOverlap(QueryRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getSequence(); // Get the forward string of read1
	UINT64 hashStringLength = this->omegaHashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement(); // Get the string from read2 according to orient.
	if(orient == 0 || orient == 2)		// orient 0
										//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
										//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
										//				OR
										// orient 2
										//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
										//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(string1.length()-start<this->minimumOverlapLength)return false;
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
//		return string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)); // If the remaining strings match.
		return strcmp(string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)).c_str(), string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)).c_str())==0; // If the remaining strings match.

	}
	else								// orient 1
										//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
										//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
										//				OR
										// orient 3
										//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
										//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(start+hashStringLength<this->minimumOverlapLength)return false;
		if(string2.length()-hashStringLength < start)
			return false;
//		return string1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
		return strcmp(string1.substr(0, start).c_str(), string2.substr(string2.length()-hashStringLength-start, start).c_str())==0; // If the remaining strings match.

	}
}

//either subject is contained in query or query is contained in subject
bool OmegaGraphConstructor::isContainedAlignment(Alignment * subjectAlignment)
{
	if(subjectAlignment->subjectStart<=0&&subjectAlignment->subjectEnd>=subjectAlignment->queryEnd)
		return true;
	else if(subjectAlignment->subjectStart>=0&&subjectAlignment->subjectEnd<=subjectAlignment->queryEnd)
		return true;
	else return false;
}
/*
bool OmegaGraphConstructor::start(bool transitiveEdgeRemoval)
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
			this->searchHashTableFullScan(subjectEdge);
			if(subjectEdge->alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				bool newalignadded=false;
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->alignmentList->size(); j++)
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


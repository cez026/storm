/*
 * QueryDatasetFilter.cpp
 *
 *  Created on: Feb 23, 2015
 *      Author: qy2
 */

#include "QueryDatasetFilter.h"

QueryDatasetFilter::QueryDatasetFilter(OmegaHashTable * omegaHashTable) {
	// TODO Auto-generated constructor stub
	this->omegaHashTable = omegaHashTable;
	this->dataset = this->omegaHashTable->getDataset();

//	this->tagList =  NULL;

//	this->superNameList = NULL;

//	this->superReadLength = NULL;
}

QueryDatasetFilter::~QueryDatasetFilter() {
	// TODO Auto-generated destructor stub
//	if(this->tagList!=NULL)
//	{
//		this->tagList->clear();
//		delete this->tagList;
//	}
//	if(this->superNameList!=NULL)
//	{
//		this->superNameList->clear();
//		delete this->superNameList;
//	}
//	if(this->superReadLength!=NULL)
//	{
//		this->superReadLength->clear();
//		delete this->superReadLength;
//	}
//	this->dataset = NULL;
}

bool QueryDatasetFilter::start()
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
			this->searchHashTable(subjectEdge); //and also mark the contained read in the query read structure
			delete subjectEdge->subjectRead;
			delete subjectEdge;
			subjectEdgesList->at(i) = NULL;
		}

	}
		// end of omp parallel


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

bool QueryDatasetFilter::start_backup()
{
	string outFileName = Config::outputfilename;
	filePointer.open(outFileName.c_str());
	CLOCKSTART;
	MEMORYSTART;
/*
	if(this->tagList==NULL)
	{
		this->tagList = new vector<UINT8>;
		this->tagList->reserve(dataset->getNumberOfUniqueReads()+1);
	}
	if(this->superNameList==NULL)
	{
		this->superNameList = new vector<string>;
		this->superNameList->reserve(dataset->getNumberOfUniqueReads()+1);
	}
	if(this->superReadLength==NULL)
	{
		this->superReadLength = new vector<int>;
		this->superReadLength->reserve(dataset->getNumberOfUniqueReads()+1);
	}


	string emptystring = "";
	for(UINT64 i = 0; i <= dataset->getNumberOfUniqueReads(); i++) // Initialization
	{
		this->tagList->push_back(0);
		this->superNameList->push_back(emptystring);
		this->superReadLength->push_back(0);
	}
*/
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
			this->searchHashTable(subjectEdge); //and also mark the contained read in the query read structure
			delete subjectEdge->subjectRead;
			delete subjectEdge;
		}

	}
		// end of omp parallel


	//clear the unused subject reads
/*		for(UINT64 i = 0; i < subjectEdgeList->size(); i++)
		{
			SubjectEdge * subjectEdge = subjectEdgeList->at(i);
			SubjectRead * thisSubjectread = subjectEdge->subjectRead;
			if(subjectEdge->DuplicateReadList!=NULL)
			{
				for(UINT16 j = 0; j < subjectEdge->DuplicateReadList->size(); j++)
				{
					 QueryRead * thisQueryRead=  subjectEdge->DuplicateReadList->at(j);
					 if(thisSubjectread->getName()>=thisQueryRead->getName())
						 {
						 	 UINT64 queryID = thisQueryRead->getIdentifier();
						 	 UINT8 tag = this->tagList->at(queryID);
						 	 this->tagList->at(queryID)= tag|0X01;
						 }
					 thisQueryRead->setFrequency(thisQueryRead->getFrequency()+1);
				}
			}


			if(subjectEdge->contained_alignmentList==NULL)
			{
				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}
			else
			{
				//add alignments to the query reads before the edges are destroyed.
				for(UINT16 j = 0; j < subjectEdge->contained_alignmentList->size(); j++)
				{
					 ContainedAlignment * alignment= subjectEdge->contained_alignmentList->at(j);

				 	 UINT64 queryID = alignment->queryRead->getIdentifier();
				 	 UINT8 tag = this->tagList->at(queryID);
				 	 this->tagList->at(queryID)= tag|0X02;
				 	 if((int)alignment->subjectRead->getReadLength() >  this->superReadLength->at(queryID))
				 	 {
				 		this->superReadLength->at(queryID)=alignment->subjectRead->getReadLength();
				 		this->superNameList->at(queryID) = alignment->subjectRead->getName();
				 	 }
				 	 else if((int)alignment->subjectRead->getReadLength() == this->superReadLength->at(queryID))
				 	 {
				 		 if(this->superNameList->at(queryID) < alignment->subjectRead->getName())
				 			this->superNameList->at(queryID) = alignment->subjectRead->getName();
				 	 }
				 	 delete alignment;

				}

				delete subjectEdge->subjectRead;
				subjectEdge->subjectRead = NULL;
			}

			delete subjectEdge;

		}
*/


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

void QueryDatasetFilter::printToFile()
{
	CLOCKSTART;
	MEMORYSTART;
	UINT64 count = 0;
	if(Config::useID == false)
	{
		for(UINT64 i = 1; i <= dataset->getNumberOfUniqueReads(); i++)
		{
			QueryRead *queryRead = dataset->getReadFromID(i);
//			int tag = this->tagList->at(i);
//			if(tag==0)
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
	//			int tag = this->tagList->at(i);
	//			if(tag==0)
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

bool QueryDatasetFilter::searchHashTable(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
	for(UINT64 j = 0; j<=subjectRead->getReadLength()-this->omegaHashTable->getHashStringLength(); j++)
	{
//		if(subjectRead->getReadLength()==35)
//		cout<<"name: "<<subjectRead->getName()<<" len "<<subjectRead->getReadLength()<<" j "<<j<<" key "<<this->omegaHashTable->getHashStringLength()<<endl;
		subString = subjectReadString.substr(j,this->omegaHashTable->getHashStringLength());
		vector<UINT64> * listOfReads=this->omegaHashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
//				UINT16 overlapOffset;
				UINT64 orientation = (data >> 62);
				QueryRead *queryRead = this->omegaHashTable->getDataset()->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.

				if(queryRead->flag4Removal==true) continue;
				if(subjectRead->getName()==queryRead->getName())continue;
				if(subjectRead->getReadLength()==queryRead->getReadLength())
				{
//					if(checkIdenticalRead(subjectRead,queryRead,orientation,j))
//					subjectEdge->addDuplicateList(queryRead);
					if(checkIdenticalRead(subjectRead,queryRead,orientation,j))
					{
						if(queryRead->getName()<subjectRead->getName())
						{
							#pragma omp critical(markContainedToQueryRead)
							 {
								 queryRead->flag4Removal = true;
							 }
						}

					}

				}
				else if(checkOverlapForContainedRead(subjectRead,queryRead,orientation,j)) // Both read need to be non contained.
				{
//YAO				int orientation = data >> 62;
//YAO					cout<<sRead->getName() <<" "<< " >>> " << qRead->getName() <<" orientation: "<<orientation<<" position: "<<j<<endl;

//YAO					cout << "S  " << sRead->getSequence() << endl;
//YAO				cout << "Q  " << qRead->getSequence() << endl;
//YAO					std::stringstream sstm;
//YAO					if(sRead->getSequence()<qRead->getSequence())
//YAO					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<sRead->getSequence()<<"||||"<<qRead->getSequence()<<endl;
//YAO					else
//YAO					sstm<<" orientation: "<<orientation<<" position: "<<j<<"||||"<<qRead->getSequence()<<"||||"<<sRead->getSequence()<<endl;
//YAO				filePointer<<sstm.str();
					if(queryRead->getReadLength()<subjectRead->getReadLength())
					{
						#pragma omp critical(markContainedToQueryRead)
						 {
							 queryRead->flag4Removal = true;
						 }
					}
//					ContainedAlignment *subjectAlignment = new ContainedAlignment(subjectRead, queryRead);
/*					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
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
					*/
//					subjectEdge->addContainedAlignment(subjectAlignment);

				}
			}
		}
	}
	return true;
}
bool QueryDatasetFilter::checkIdenticalRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	//containing redundancy
//	if(start!=0 || start!=read1->getReadLength()-this->omegaHashTable->getHashStringLength())return false;
	if(start!=0)return false;
	string string1=read1->getSequence();
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement();
//	if(string1 == string2) return true;
	if(strcmp(string1.c_str(), string2.c_str())==0) return true;
	else return false;
}



/**********************************************************************************************************************
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	orient 0 means prefix of forward of the read2
	orient 1 means suffix of forward of the read2
	orient 2 means prefix of reverse of the read2
	orient 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read2 is contained in read1.
**********************************************************************************************************************/
bool QueryDatasetFilter::checkOverlapForContainedRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getSequence(); // Get the forward of read1
	UINT64 hashStringLength = this->omegaHashTable->getHashStringLength(), lengthRemaining1, lengthRemaining2;
	string string2 = (orient == 0 || orient== 1) ? read2->getSequence() : read2->reverseComplement(); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)
									// orient 0
									//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
									//				OR
									// orient 2
									//	 >---*****MMMMMMMMMMMMMMM*******------> read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	{
		lengthRemaining1 = string1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
//			cout<<"UP: "<<"read1 "<<read1->getName()<<" len "<<read1->getReadLength()<<"substr: "<<start + hashStringLength<<"||"<<lengthRemaining2<<endl;
//			cout<<"UP: "<<"read2 "<<read2->getName()<<" len "<<read2->getReadLength()<<"substr: "<<hashStringLength<<"||"<<lengthRemaining2<<endl;
//			cout<<"---"<<string1.substr(start + hashStringLength, lengthRemaining2)<<endl;
//			cout<<"---"<<string2.substr(hashStringLength, lengthRemaining2)<<endl;

//			return string1.substr(start + hashStringLength, lengthRemaining2) == string2.substr(hashStringLength, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
			return strcmp(string1.substr(start + hashStringLength, lengthRemaining2).c_str(), string2.substr(hashStringLength, lengthRemaining2).c_str())==0; // If the remaining of the string match, then read2 is contained in read1

		}
	}
	else							// orient 1
									//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
									//				OR
									// orient 3
									//	 >---*****MMMMMMMMMMMMMMM-------------> read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
//			cout<<"DOWN: "<<"read1 "<<read1->getName()<<" len "<<read1->getReadLength()<<"substr: "<<start - lengthRemaining2<<"||"<<lengthRemaining2<<endl;
//			cout<<"DOWN: "<<"read2 "<<read2->getName()<<" len "<<read2->getReadLength()<<"substr: "<<'0'<<"||"<<lengthRemaining2<<endl;

			return strcmp(string1.substr(start - lengthRemaining2, lengthRemaining2).c_str(), string2.substr(0, lengthRemaining2).c_str())==0; // If the remaining of the string match, then read2 is contained in read1
		}
	}
	return false;

}

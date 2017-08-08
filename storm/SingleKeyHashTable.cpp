/*
 * SingleKeyHashTable.cpp
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#include "SingleKeyHashTable.h"

SingleKeyHashTable::SingleKeyHashTable() {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength = Config::hashKeyLength;
	this->maxMismatch = Config::maxMismatch;
	this->maxErrorRate = ((double)(Config::maxErrorRate))/100;
	this->numberOfMode = 4;
	this->numberOfMSB = 2;
	this->numberOfLSB = 62;
	this->dataSet = NULL;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

SingleKeyHashTable::SingleKeyHashTable(QueryDataset * qDataset) {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength = Config::hashKeyLength;
	this->maxMismatch = Config::maxMismatch;
	this->maxErrorRate = ((double)(Config::maxErrorRate))/100;
	this->numberOfMode = 4;
	this->numberOfMSB = 2;
	this->numberOfLSB = 62;
	this->dataSet = qDataset;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

SingleKeyHashTable::~SingleKeyHashTable() {
	// TODO Auto-generated destructor stub
	if(this->hashTable!=NULL)
	delete hashTable;
	this->dataSet = NULL;//we don't delete dataSet here
}

bool SingleKeyHashTable::createHashTables()
{
	if(this->dataSet==NULL)
	{
		cout<<"no data set"<<endl;
		return false;
	}
	else if(this->hashTable!=NULL)
	{
		cout<<"Hash Table already exists."<<endl;
		return false;
	}
	else
	{
		this->hashTable = new HashTable(8*this->dataSet->getNumberOfUniqueReads());

		return true;
	}
}

// 00 = 0 means prefix of the forward string.
// 01 = 1 means suffix of the forward string.
// 10 = 2 means prefix of the reverse string.
// 11 = 3 means suffix of the reverse string.
string SingleKeyHashTable::getReadSubstring(UINT64 readID, UINT8 mode)
{
	QueryRead * read = this->dataSet->getReadFromID(readID);
	string str = (mode == 0 || mode == 1) ? read->getSequence() : read->reverseComplement();
	string subStr = (mode == 0 || mode == 2) ? str.substr(0,this->hashKeyLength) : str.substr(str.length() - this->hashKeyLength, this->hashKeyLength);

	return subStr;
}


bool SingleKeyHashTable::insertQueryDataset(QueryDataset* d)
{
	if(this->hashTable==NULL)
	{
		cout<<"Hash Table hasn't been created yet."<<endl;
		return false;
	}
	else
	{
		UINT64 datasetsize = this->dataSet->getNumberOfUniqueReads();
		UINT64 currentID = 1;
		while(currentID<=datasetsize)
		{
			if(currentID%1000000 == 0)
				cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
			QueryRead * read = this->dataSet->getReadFromID(currentID);
			string forwardRead = read->getSequence();
			string reverseRead = read->reverseComplement();
			string prefixForward = forwardRead.substr(0,this->hashKeyLength); 											// Prefix of the forward string.
			string suffixForward = forwardRead.substr(forwardRead.length() - this->hashKeyLength,this->hashKeyLength);	// Suffix of the forward string.
			string prefixReverse = reverseRead.substr(0,this->hashKeyLength);											// Prefix of the reverse string.
			string suffixReverse = reverseRead.substr(reverseRead.length() - this->hashKeyLength,this->hashKeyLength);	// Suffix of the reverse string.
			insertQueryRead(read, prefixForward, 0);
			insertQueryRead(read, suffixForward, 1);
			insertQueryRead(read, prefixReverse, 2);
			insertQueryRead(read, suffixReverse, 3);
			currentID++;
		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;



		return true;

	}
}

bool SingleKeyHashTable::insertQueryDataset_leftReadForMerge(QueryDataset* d)
{
	if(this->hashTable==NULL)
	{
		cout<<"Hash Table hasn't been created yet."<<endl;
		return false;
	}
	else
	{
		UINT64 datasetsize = this->dataSet->getNumberOfReads();
		UINT64 currentID = 1;
		while(currentID<=datasetsize)
		{
			if(currentID%1000000 == 0)
				cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
			QueryRead * read = this->dataSet->getReadFromID(currentID);
			if(read->getPairPosition()==1)
			{
			string forwardRead = read->getSequence();
			string reverseRead = read->reverseComplement();
			string suffixForward = forwardRead.substr(forwardRead.length() - this->hashKeyLength,this->hashKeyLength);	// Suffix of the forward string.
			string prefixReverse = reverseRead.substr(0,this->hashKeyLength);											// Prefix of the reverse string.
			insertQueryRead(read, suffixForward, 1);
			insertQueryRead(read, prefixReverse, 2);
			}
			currentID++;

		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;



		return true;

	}
}

bool SingleKeyHashTable::insertQueryDataset_rightAlignOnly(QueryDataset* d)
{
	if(this->hashTable==NULL)
	{
		cout<<"Hash Table hasn't been created yet."<<endl;
		return false;
	}
	else
	{
		UINT64 datasetsize = this->dataSet->getNumberOfReads();
		UINT64 currentID = 1;
		while(currentID<=datasetsize)
		{
			if(currentID%1000000 == 0)
				cout << setw(10) << currentID << " reads inserted in the hash table. " << endl;
			QueryRead * read = this->dataSet->getReadFromID(currentID);

			string forwardRead = read->getSequence();
			string reverseRead = read->reverseComplement();
			string suffixForward = forwardRead.substr(forwardRead.length() - this->hashKeyLength,this->hashKeyLength);	// Suffix of the forward string.
			string prefixReverse = reverseRead.substr(0,this->hashKeyLength);											// Prefix of the reverse string.
			insertQueryRead(read, suffixForward, 1);
			insertQueryRead(read, prefixReverse, 2);
			currentID++;

		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;



		return true;

	}
}

bool SingleKeyHashTable::insertQueryRead(QueryRead *read, string subString, UINT8 mode)
{
	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);
	while(!hashTable->isEmptyAt(index))
	{
		vector<UINT64>* readList = this->hashTable->getReadIDListAt(index);
		UINT64 data = readList->at(0);
		UINT64 keyreadID = data & 0X3FFFFFFFFFFFFFFF;
		UINT64 keymode = data >> this->numberOfLSB;


		string keyStr = this->getReadSubstring(keyreadID,keymode);


//		if(keyStr == subString)
//				break;
		if(strcmp(keyStr.c_str(),subString.c_str())==0)
				break;
		numberOfHashCollision++;
		currentCollision++;
		index = (index == hashTable->getHashTableSize() - 1) ? 0: index + 1; 	// Increment the index
	}
	UINT64 a = ((UINT64)(mode) << this->numberOfLSB);
	UINT64 b = read->getIdentifier();
	UINT64 combinedID = a | b;
	hashTable->insertReadIDAt(index,combinedID);							// Add the string in the list.

	if(currentCollision> this->maxSingleHashCollision)
		this->maxSingleHashCollision = currentCollision;
//	if(currentCollision > 1000)
//	{
//		cout << currentCollision << " collisions for read " << read->getIdentifier() << " " << subString << endl;
//	}
	return true;
}

vector<UINT64> * SingleKeyHashTable::getListOfReads(string subString)
{

	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);	// Get the index using the hash function.
	while(!hashTable->isEmptyAt(index))
	{
		vector<UINT64>* readList = this->hashTable->getReadIDListAt(index);
		UINT64 data = readList->at(0);
		UINT64 keyreadID = data & 0X3FFFFFFFFFFFFFFF;
		UINT64 keymode = data >> this->numberOfLSB;
		string keyStr = this->getReadSubstring(keyreadID,keymode);


//		if(keyStr == subString)
//				break;
		if(strcmp(keyStr.c_str(),subString.c_str())==0)
				break;
		currentCollision++;
		if(currentCollision>this->maxSingleHashCollision)return NULL;
		index = (index == hashTable->getHashTableSize() - 1) ? 0: index + 1; 	// Increment the index
	}

	return hashTable->getReadIDListAt(index);	// return the index.
}



bool SingleKeyHashTable::searchHashTable(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength; j++)
	{
		subString = subjectReadString.substr(j,this->hashKeyLength);
		vector<UINT64> * listOfReads=this->getListOfReads(subString); // Search the string in the hash table.
		if(listOfReads!=NULL && !listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(INT64 k = 0; k < (INT64)listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT64 queryReadID = data & 0X3FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID); 	// Least significant 62 bits store the read number.
				if(subjectRead->getName()<=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
				{
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}

				if(subjectRead->flag4Removal==false && queryRead->flag4Removal==false) // Both read need to be non contained.
				{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);

					if(this->createAlignment(subjectAlignment, queryMode, j))
					{

					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
			}
		}
	}
	return true;
}

bool SingleKeyHashTable::searchHashTableFullScan(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength; j++)
	{
		subString = subjectReadString.substr(j,this->hashKeyLength);
		vector<UINT64> * listOfReads=this->getListOfReads(subString); // Search the string in the hash table.
		if(listOfReads!=NULL && !listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(INT64 k = 0; k < (INT64)listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT64 queryReadID = data & 0X3FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID); 	// Least significant 62 bits store the read number.
				if(subjectRead->getName()==queryRead->getName()) 			// This is different from half scan.
				{
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}

				if(subjectRead->flag4Removal==false && queryRead->flag4Removal==false) // Both read need to be non contained.
				{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);

					if(this->createAlignment(subjectAlignment, queryMode, j))
					{

					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
			}
		}
	}
	return true;
}

bool SingleKeyHashTable::searchHashTable(SubjectEdge * subjectEdge, bool fullScan)
{
	if(fullScan == false) return this->searchHashTable(subjectEdge);
	else return this->searchHashTableFullScan(subjectEdge);
}


bool SingleKeyHashTable::searchHashTableFullScanAndMarkQueryRead(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength; j++)
	{
		subString = subjectReadString.substr(j,this->hashKeyLength);
		vector<UINT64> * listOfReads=this->getListOfReads(subString); // Search the string in the hash table.
		if(listOfReads!=NULL && !listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(INT64 k = 0; k < (INT64)listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT64 queryReadID = data & 0X3FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID); 	// Least significant 62 bits store the read number.
				if(subjectRead->getName()==queryRead->getName()) 			// This is different from half scan.
				{
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}

				if(queryRead->flag4Removal==false) // Both read need to be non contained.
				{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);

					if(this->createAlignment(subjectAlignment, queryMode, j))
					{
						if(subjectAlignment->subjectStart<=0&&subjectAlignment->subjectEnd>=subjectAlignment->queryEnd)
						{
							if(subjectRead->getReadLength()==queryRead->getReadLength())//identical reads
							{
								if(queryRead->getName()<subjectRead->getName()) //only keep the larger name
								{
//									queryRead->setRemovalTag(true);
//									#pragma omp critical(markContainedToQueryRead)
//									 {
										 queryRead->flag4Removal = true;
//									 }
								}
							}
							else //contained reads
							{
//								queryRead->setRemovalTag(true);
//								#pragma omp critical(markContainedToQueryRead)
//								 {
									 queryRead->flag4Removal = true;
//								 }
							}
						}

					}

					//no alignment needs to be stored
					delete subjectAlignment;
				}
			}
		}
	}
	return true;
}

vector<UINT64> * SingleKeyHashTable::getListOfReads(string subString, UINT64& hashValue)
{

	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);	// Get the index using the hash function.
	hashValue = index;
	while(!hashTable->isEmptyAt(index))
	{
		vector<UINT64>* readList = this->hashTable->getReadIDListAt(index);
		UINT64 data = readList->at(0);
		UINT64 keyreadID = data & 0X3FFFFFFFFFFFFFFF;
		UINT64 keymode = data >> this->numberOfLSB;
		string keyStr = this->getReadSubstring(keyreadID,keymode);


//		if(keyStr == subString)
//				break;
		if(strcmp(keyStr.c_str(),subString.c_str())==0)
				break;
		currentCollision++;
		if(currentCollision>this->maxSingleHashCollision)return NULL;
		index = (index == hashTable->getHashTableSize() - 1) ? 0: index + 1; 	// Increment the index
	}

	return hashTable->getReadIDListAt(index);	// return the index.
}

bool SingleKeyHashTable::searchHashTableFullScanAndRecordKeys(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength; j++)
	{
		subString = subjectReadString.substr(j,this->hashKeyLength);
		UINT64 hashValue;
		vector<UINT64> * listOfReads=this->getListOfReads(subString, hashValue); // Search the string in the hash table.
		subjectEdge->addKeyToMap(hashValue,j);
		if(listOfReads!=NULL && !listOfReads->empty()) // If there are some reads that contain subString as prefix or suffix of the read or their reverse complement
		{
			for(INT64 k = 0; k < (INT64)listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT64 queryReadID = data & 0X3FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID); 	// Least significant 62 bits store the read number.
				if(subjectRead->getName()==queryRead->getName()) 			// This is different from half scan.
				{
					continue;
				}
				else
				{
//					cout << sRead->getName() << " >>> " << qRead->getName() << " ????? " <<j<<" and " << k << endl;

				}

				if(subjectRead->flag4Removal==false && queryRead->flag4Removal==false) // Both read need to be non contained.
				{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);

					if(this->createAlignment(subjectAlignment, queryMode, j))
					{

					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
			}
		}
	}
	return true;
}

bool SingleKeyHashTable::createAlignment(Alignment* subjectAlignment, UINT8 queryMode, UINT64 subjectKeyStart)
{
	QueryRead* queryRead = subjectAlignment->queryRead;
	SubjectRead* subjectRead = subjectAlignment->subjectRead;
	string subjectString = subjectRead->getSequence();
	switch (queryMode) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
	{
	case 0:
	{
		subjectAlignment->queryOrientation = true;
		subjectAlignment->subjectStart = -subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		string queryString = queryRead->getSequence();
		int remainStart = this->hashKeyLength;
		int remainEnd = (subjectAlignment->subjectEnd <= subjectAlignment->queryEnd)?subjectAlignment->subjectEnd:subjectAlignment->queryEnd;
		bool alignsucess = this->doAlignment(subjectAlignment, queryString, subjectString, remainStart, remainEnd);
		return alignsucess;
		break;
	}
	case 1:
	{
		subjectAlignment->queryOrientation = true;
		subjectAlignment->subjectStart = queryRead->getReadLength()-this->hashKeyLength-subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		string queryString = queryRead->getSequence();
		int remainStart = (0>=subjectAlignment->subjectStart)?0:subjectAlignment->subjectStart;
		int remainEnd = subjectAlignment->queryEnd - this->hashKeyLength;
		bool alignsucess = this->doAlignment(subjectAlignment, queryString, subjectString, remainStart, remainEnd);
		return alignsucess;
		break;
	}
	case 2:
	{
		subjectAlignment->queryOrientation = false;
		subjectAlignment->subjectStart = -subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		string queryString = queryRead->reverseComplement();
		int remainStart = this->hashKeyLength;
		int remainEnd = (subjectAlignment->subjectEnd <= subjectAlignment->queryEnd)?subjectAlignment->subjectEnd:subjectAlignment->queryEnd;
		bool alignsucess = this->doAlignment(subjectAlignment, queryString, subjectString, remainStart, remainEnd);
		return alignsucess;
		break;
	}
	case 3:
	{
		subjectAlignment->queryOrientation = false;
		subjectAlignment->subjectStart = queryRead->getReadLength()-this->hashKeyLength-subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		string queryString = queryRead->reverseComplement();
		int remainStart = (0>=subjectAlignment->subjectStart)?0:subjectAlignment->subjectStart;
		int remainEnd = subjectAlignment->queryEnd - this->hashKeyLength;
		bool alignsucess = this->doAlignment(subjectAlignment, queryString, subjectString, remainStart, remainEnd);
		return alignsucess;
		break;
	}
	default: return false;
	}

}



bool SingleKeyHashTable::doAlignment(Alignment* subjectAlignment,string& queryString, string& subjectString, int remainStart, int remainEnd)
{
	int overlaplength = remainEnd-remainStart+1+this->hashKeyLength;
	if(overlaplength>=this->minimumOverlapLength)
	{
		int maxError;
		if(Config::useErrorRate==true)
		{
			maxError = floor(overlaplength*this->maxErrorRate);

		}
		else maxError = this->maxMismatch;

		int currentMismatchCount = 0;
		for(int i=remainStart;i<=remainEnd;i++)
		{
			char queryBase = queryString.at(i);
			char subjectBase = subjectString.at(i-subjectAlignment->subjectStart);
			if(queryBase!=subjectBase)
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxError)return false;
				subjectAlignment->insertSubstitution(i, subjectBase);
			}
		}
		return true;
	}
	else return false;
}


/*

bool SingleKeyHashTable::doAlignment(Alignment* subjectAlignment,string& queryString, string& subjectString, int remainStart, int remainEnd)
{
	bool successflag = true;
	int overlaplength = remainEnd-remainStart+1+this->hashKeyLength;
	if(overlaplength>=this->minimumOverlapLength)
	{
		int maxError;
		if(Config::useErrorRate==true)
		{
			maxError = floor(overlaplength*this->maxErrorRate);

		}
		else maxError = this->maxMismatch;

		int currentMismatchCount = 0;
		int cursor = remainStart;
		while(cursor<=remainEnd&&successflag==true)
		{
			if(currentMismatchCount==maxError)
			{
				successflag = (queryString.substr(cursor,remainEnd-cursor)==subjectString.substr(cursor-subjectAlignment->subjectStart,remainEnd-cursor));
				cursor=remainEnd+1;
			}
			else if(currentMismatchCount<maxError)
			{
				char queryBase = queryString.at(cursor);
				char subjectBase = subjectString.at(cursor-subjectAlignment->subjectStart);
				if(queryBase!=subjectBase)
				{
					currentMismatchCount++;
					subjectAlignment->insertSubstitution(cursor, subjectBase);
				}
				cursor++;
			}
			else
			{
				successflag = false;
			}
		}

	}
	else successflag = false;

	return successflag;
}


bool SingleKeyHashTable::doAlignment(Alignment* align, string mode, int subjectStart)
{

	if(mode=="forwardprefix")
	{
		align->queryOrientation = true;
		if(align->subjectReadSequence.length()- subjectStart - hashKeyLength >= align->queryRead->getSequence().length() - hashKeyLength) // The overlap must continue till the end.
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(subjectStart + hashKeyLength, align->subjectReadSequence.length()-(subjectStart + hashKeyLength));
		string restQuery = align->queryRead->getSequence().substr(hashKeyLength,  restSubject.length());
		int currentMismatchCount=0;
		for(unsigned int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
			}
		}
		align->subjectStart = -subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "forwardsuffix")
	{
		align->queryOrientation = true;
		if(align->queryRead->getSequence().length()-hashKeyLength <= subjectStart)
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(0, subjectStart);
		string restQuery = align->queryRead->getSequence().substr(align->queryRead->getReadLength()-hashKeyLength-subjectStart, restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(align->queryRead->getReadLength()-hashKeyLength-subjectStart+i, restSubject.at(i)));
			}
		}
		align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "reverseprefix")
	{
		align->queryOrientation = false;
		if(align->subjectReadSequence.length()- subjectStart - hashKeyLength >= align->queryRead->getSequence().length() - hashKeyLength) // The overlap must continue till the end.
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(subjectStart + hashKeyLength, align->subjectReadSequence.length()-(subjectStart + hashKeyLength));
		string restQuery = align->queryRead->reverseComplement().substr(hashKeyLength,  restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
			}
		}
		align->subjectStart = -subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else if(mode == "reversesuffix")
	{
		align->queryOrientation = false;
		if(align->queryRead->getSequence().length()-hashKeyLength <= subjectStart)
			return checkForContainedAlignment(align, mode, subjectStart);
		string restSubject = align->subjectReadSequence.substr(0, subjectStart);
		string restQuery = align->queryRead->reverseComplement().substr(align->queryRead->getReadLength()-hashKeyLength-subjectStart, restSubject.length());
		int currentMismatchCount=0;
		for(int i=0; i<restQuery.length();i++)
		{
			if(restQuery.at(i)!=restSubject.at(i))
			{
				currentMismatchCount++;
				if(currentMismatchCount>maxMismatch)return false;
				align->editInfor.insert(std::pair<int, char>(align->queryRead->getReadLength()-hashKeyLength-subjectStart+i, restSubject.at(i)));
			}
		}
		align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
		align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
		align->queryEnd = align->queryRead->getReadLength()-1;
	}
	else return false;

	return true;
}

//the choice of the start and stop position should meet the minimum overlap requirement.
bool SingleKeyHashTable::subjectWindowRange(int& startpoint, int& stoppoint, string mode, string& subjectRead)
{
	if(mode=="forwardprefix" || mode == "reverseprefix")
	{
		startpoint = 0;
		stoppoint = subjectRead.length()-this->minimumOverlapLength;
	}
	else if(mode == "forwardsuffix" || mode == "reversesuffix")
	{
		startpoint = this->minimumOverlapLength - hashKeyLength;
		stoppoint = subjectRead.length()-hashKeyLength;

	}
	else return false;

	return true;
}

bool SingleKeyHashTable::checkForContainedAlignment(Alignment* align, string mode, int subjectStart)
{
	string subjectString=align->subjectReadSequence; // Get the forward of read1
	string queryString="";
//	string queryString = (mode=="forwardprefix" || mode=="forwardsuffix") ? align->queryRead->getSequence() : align->queryRead->reverseComplement(); // Get the string in read2 based on the orientation.
	if(mode=="forwardprefix" || mode=="forwardsuffix")
	{
		queryString = align->queryRead->getSequence();
		align->queryOrientation = true;
	}
	else if(mode=="reverseprefix" || mode=="reversesuffix")
	{
		queryString = align->queryRead->reverseComplement();
		align->queryOrientation = false;
	}
	else return false;

	if(mode=="forwardprefix" || mode=="reverseprefix")
									// mode = forwardprefix
									//   >--------MMMMMMMMMMMMMMM*******------> subject read1      M means match found by hash table
									//            MMMMMMMMMMMMMMM*******>       query read2      * means we need to check these characters for match
									//				OR
									// mode = reverseprefix
									//	 >---*****MMMMMMMMMMMMMMM*******------> subject read1
									//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of query read2
	{
		int restSubjectLength = subjectString.length() - subjectStart - this->hashKeyLength; 	// This is the remaining of read1
		int restQueryLength = queryString.length() - this->hashKeyLength; 	// This is the remaining of read2
		if(restSubjectLength >= restQueryLength)
		{
			string restSubject = subjectString.substr(subjectStart + hashKeyLength, restQueryLength);
			string restQuery = queryString.substr(hashKeyLength,  restQueryLength);
			int currentMismatchCount=0;
			for(int i=0; i<restQuery.length();i++)
			{
				if(restQuery.at(i)!=restSubject.at(i))
				{
					currentMismatchCount++;
					if(currentMismatchCount>maxMismatch)return false;
					align->editInfor.insert(std::pair<int, char>(hashKeyLength+i, restSubject.at(i)));
				}
			}
			align->subjectStart = -subjectStart;
			align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
			align->queryEnd = align->queryRead->getReadLength()-1;

		}
	}
	else if(mode=="forwardsuffix" || mode=="reversesuffix")
									// mode = forwardsuffix
									//   >---*****MMMMMMMMMMMMMMM-------------> subject read1      M means match found by hash table
									//      >*****MMMMMMMMMMMMMMM       		query read2      * means we need to check these characters for match
									//				OR
									// mode = reversesuffix
									//	 >---*****MMMMMMMMMMMMMMM-------------> subject read1
									//		<*****MMMMMMMMMMMMMMM				Reverse Complement of query Read2
	{
		int restSubjectLength = subjectStart;
		int restQueryLength = queryString.length() - this->hashKeyLength;
		if(restSubjectLength >= restQueryLength)
		{
			string restSubject = subjectString.substr(subjectStart-restQueryLength, restQueryLength);
			string restQuery = queryString.substr(0, restQueryLength);
			int currentMismatchCount=0;
			for(int i=0; i<restQuery.length();i++)
			{
				if(restQuery.at(i)!=restSubject.at(i))
				{
					currentMismatchCount++;
					if(currentMismatchCount>maxMismatch)return false;
					align->editInfor.insert(std::pair<int, char>(i, restSubject.at(i)));
				}
			}
			align->subjectStart = align->queryRead->getReadLength()-hashKeyLength-subjectStart;
			align->subjectEnd = align->subjectStart + align->subjectReadSequence.length()-1;
			align->queryEnd = align->queryRead->getReadLength()-1;
		}
	}
	else return false;

	return true;

}

bool SingleKeyHashTable::singleKeySearch(edge & Edge)
{

	string subjectRead = Edge.subjectReadSequence;


	for(int i =0; i<this->hashTableNameList.size();i++)
	{


	string modestring = this->hashTableNameList.at(i);

	int startpoint,stoppoint;
	if(!subjectWindowRange(startpoint, stoppoint, modestring, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subString = subjectRead.substr(j, hashKeyLength);
	vector<UINT64> * currentIDList = hashTableMap.at(modestring)->getReadIDListOfReads(subString);
	for(UINT64 k=0;currentIDList!=NULL&&k<currentIDList->size();k++)
	{
		UINT64 currentID = currentIDList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName>=Edge.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		//however, this will cause missing of the contained read detection/alignment. Accept this reality.
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = Edge.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(doAlignment(align, modestring, j)==false) delete align;
			else Edge.alignmentList.push_back(align);

		}
	}
	}

	}

}

bool SingleKeyHashTable::singleKeySearch(SubjectAlignment & subjectAlign)
{

	string subjectRead = subjectAlign.subjectReadSequence;


	for(int i =0; i<this->hashTableNameList.size();i++)
	{


	string modestring = this->hashTableNameList.at(i);

	int startpoint,stoppoint;
	if(!subjectWindowRange(startpoint, stoppoint, modestring, subjectRead)) return false;//guarantee it meets the minimum overlaplength requirement for alignments.

	for(int j=startpoint;j<=stoppoint;j++)
	{
	string subString = subjectRead.substr(j, hashKeyLength);
	vector<UINT64>* currentIDList = hashTableMap.at(modestring)->getReadIDListOfReads(subString);
	for(UINT64 k=0;currentIDList!=NULL&&k<currentIDList->size();k++)
	{
		UINT64 currentID = currentIDList->at(k);
		QueryRead* queryRead = queryDataSet->getReadFromID(currentID);
		string querySequence = queryRead->getSequence();
		string queryName = queryRead->getName();
		bool alignFlag = true;
		if(queryName==subjectAlign.subjectReadName)alignFlag = false; //only align when queryName is alphabetical smaller than subject, removing half of the pair-wise alignment redundancy
		//however, this will cause missing of the contained read detection/alignment. Accept this reality.
		if(alignFlag)
		{
			Alignment* align = new Alignment();
			align->subjectReadName = subjectAlign.subjectReadName;
			align->subjectReadSequence = subjectRead;
			align->queryRead = queryRead;
			if(doAlignment(align, modestring, j)==false) delete align;
			else subjectAlign.queryAlignmentList.push_back(align);

		}
	}
	}

	}

}
*/

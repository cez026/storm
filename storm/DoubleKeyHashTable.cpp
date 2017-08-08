/*
 * DoubleKeyHashTable.cpp
 *
 *  Created on: Feb 15, 2015
 *      Author: qy2
 */

#include "DoubleKeyHashTable.h"

DoubleKeyHashTable::DoubleKeyHashTable() {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength_left = Config::hashKeyLength_left;
	this->hashKeyLength_right = Config::hashKeyLength_right;
	this->maxMismatch = Config::maxMismatch;
	this->maxErrorRate = ((double)(Config::maxErrorRate))/100;
	this->numberOfMode = 8;
	this->numberOfMSB = 3;
	this->numberOfLSB = 61;
	this->dataSet = NULL;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

DoubleKeyHashTable::DoubleKeyHashTable(QueryDataset * qDataset) {
	// TODO Auto-generated constructor stub
	this->minimumOverlapLength = Config::minimumOverlapLength;
	this->hashKeyLength_left = Config::hashKeyLength_left;
	this->hashKeyLength_right = Config::hashKeyLength_right;
	this->maxMismatch = Config::maxMismatch;
	this->maxErrorRate = ((double)(Config::maxErrorRate))/100;
	this->numberOfMode = 8;
	this->numberOfMSB = 3;
	this->numberOfLSB = 61;
	this->dataSet = qDataset;
	this->hashTable = NULL;
	numberOfHashCollision = 0;
	maxSingleHashCollision = 0;
}

DoubleKeyHashTable::~DoubleKeyHashTable() {
	// TODO Auto-generated destructor stub
	if(this->hashTable!=NULL)
	delete hashTable;
	this->dataSet = NULL;//we don't delete dataSet here
}

bool DoubleKeyHashTable::createHashTables()
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
		this->hashTable = new HashTable(16*this->dataSet->getNumberOfUniqueReads());

		return true;
	}
}

//left key
// 000 = 0 means prefix of the forward string.
// 001 = 1 means suffix of the forward string.
// 010 = 2 means prefix of the reverse string.
// 011 = 3 means suffix of the reverse string.
//right key
// 100 = 4 means prefix of the forward string.
// 101 = 5 means suffix of the forward string.
// 110 = 6 means prefix of the reverse string.
// 111 = 7 means suffix of the reverse string.
string DoubleKeyHashTable::getReadSubstring(UINT64 readID, UINT8 mode)
{
	QueryRead * read = this->dataSet->getReadFromID(readID);
	string str = (mode == 0 || mode == 1 || mode == 4 || mode == 5) ? read->getSequence() : read->reverseComplement();
	string subStr;
	switch(mode)
	{
	case 0:
	case 2:
		subStr = str.substr(0,this->hashKeyLength_left);
		break;
	case 1:
	case 3:
		subStr = str.substr(str.length() - this->hashKeyLength_left - this->hashKeyLength_right, this->hashKeyLength_left);
		break;
	case 4:
	case 6:
		subStr = str.substr(this->hashKeyLength_left,this->hashKeyLength_right);
		break;
	case 5:
	case 7:
		subStr = str.substr(str.length() - this->hashKeyLength_right, this->hashKeyLength_right);
		break;
	default:
		cout<<"wrong mode number"<<endl;
		subStr="";
		break;
	}
	return subStr;
}

bool DoubleKeyHashTable::insertQueryDataset(QueryDataset* d)
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
			string prefixForward_left = forwardRead.substr(0,this->hashKeyLength_left);
			string prefixForward_right = forwardRead.substr(this->hashKeyLength_left,this->hashKeyLength_right);
			string suffixForward_left = forwardRead.substr(forwardRead.length() - this->hashKeyLength_left - this->hashKeyLength_right,this->hashKeyLength_left);
			string suffixForward_right = forwardRead.substr(forwardRead.length() - this->hashKeyLength_right, this->hashKeyLength_right);
			string prefixReverse_left = reverseRead.substr(0,this->hashKeyLength_left);
			string prefixReverse_right = reverseRead.substr(this->hashKeyLength_left,this->hashKeyLength_right);
			string suffixReverse_left = reverseRead.substr(reverseRead.length() - this->hashKeyLength_left - this->hashKeyLength_right,this->hashKeyLength_left);
			string suffixReverse_right = reverseRead.substr(reverseRead.length() - this->hashKeyLength_right, this->hashKeyLength_right);

			insertQueryRead(read, prefixForward_left, 0);
			insertQueryRead(read, suffixForward_left, 1);
			insertQueryRead(read, prefixReverse_left, 2);
			insertQueryRead(read, suffixReverse_left, 3);
			insertQueryRead(read, prefixForward_right, 4);
			insertQueryRead(read, suffixForward_right, 5);
			insertQueryRead(read, prefixReverse_right, 6);
			insertQueryRead(read, suffixReverse_right, 7);
			currentID++;
		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;


		return true;

	}
}


bool DoubleKeyHashTable::insertQueryDataset_leftReadForMerge(QueryDataset* d)
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
			string suffixForward_left = forwardRead.substr(forwardRead.length() - this->hashKeyLength_left - this->hashKeyLength_right,this->hashKeyLength_left);
			string suffixForward_right = forwardRead.substr(forwardRead.length() - this->hashKeyLength_right, this->hashKeyLength_right);
			string prefixReverse_left = reverseRead.substr(0,this->hashKeyLength_left);
			string prefixReverse_right = reverseRead.substr(this->hashKeyLength_left,this->hashKeyLength_right);

			insertQueryRead(read, suffixForward_left, 1);
			insertQueryRead(read, prefixReverse_left, 2);
			insertQueryRead(read, suffixForward_right, 5);
			insertQueryRead(read, prefixReverse_right, 6);
			}
			currentID++;
		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;


		return true;

	}
}

bool DoubleKeyHashTable::insertQueryDataset_rightAlignOnly(QueryDataset* d)
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
			string suffixForward_left = forwardRead.substr(forwardRead.length() - this->hashKeyLength_left - this->hashKeyLength_right,this->hashKeyLength_left);
			string suffixForward_right = forwardRead.substr(forwardRead.length() - this->hashKeyLength_right, this->hashKeyLength_right);
			string prefixReverse_left = reverseRead.substr(0,this->hashKeyLength_left);
			string prefixReverse_right = reverseRead.substr(this->hashKeyLength_left,this->hashKeyLength_right);

			insertQueryRead(read, suffixForward_left, 1);
			insertQueryRead(read, prefixReverse_left, 2);
			insertQueryRead(read, suffixForward_right, 5);
			insertQueryRead(read, prefixReverse_right, 6);

			currentID++;
		}
		cout<<"Hash Table "<<" maximum collision number is: "<< this->numberOfHashCollision<<endl;
		cout<<"Hash Table "<<" maximum single read collision number is: "<< this->maxSingleHashCollision<<endl;


		return true;

	}
}


bool DoubleKeyHashTable::insertQueryRead(QueryRead *read, string subString, UINT8 mode)
{
	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);
	while(!hashTable->isEmptyAt(index))
	{
		vector<UINT64>* readList = this->hashTable->getReadIDListAt(index);
		UINT64 data = readList->at(0);
		UINT64 keyreadID = data & 0X1FFFFFFFFFFFFFFF;
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

	if(currentCollision > this->maxSingleHashCollision)
		this->maxSingleHashCollision = currentCollision;
//	if(currentCollision > 1000)
//	{
//		cout << currentCollision << " collisions for read " << read->getIdentifier() << " " << subString << endl;
//	}
	return true;
}

vector<UINT64> * DoubleKeyHashTable::getListOfReads(string subString)
{

	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);	// Get the index using the hash function.
	while(!hashTable->isEmptyAt(index))
	{
		vector<UINT64>* readList = this->hashTable->getReadIDListAt(index);
		UINT64 data = readList->at(0);
		UINT64 keyreadID = data & 0X1FFFFFFFFFFFFFFF;
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

vector<UINT64> * DoubleKeyHashTable::getListOfReads(string subString, UINT64& hashValue)
{

	UINT64 currentCollision =0;

	UINT64 index = this->hashTable->hashFunction(subString);	// Get the index using the hash function.
	hashValue = index;
	while(!hashTable->isEmptyAt(index))
	{
		vector<UINT64>* readList = this->hashTable->getReadIDListAt(index);
		UINT64 data = readList->at(0);
		UINT64 keyreadID = data & 0X1FFFFFFFFFFFFFFF;
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

bool DoubleKeyHashTable::doAlignmentWithSeed(Alignment* subjectAlignment,string& queryString, string& subjectString, int start, int end, int seedStart, int seedEnd)
{
	int overlaplength = end-start+1;
	if(overlaplength>=this->minimumOverlapLength)
	{
		int maxError;
		if(Config::useErrorRate==true)
		{
			maxError = floor(overlaplength*this->maxErrorRate);

		}
		else maxError = this->maxMismatch;
		int currentMismatchCount=0;
		int i = start;
		while(i<=end)
		{
			if(i==seedStart)
			{
				i = seedEnd;

			}
			else
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
			i++;
		}
		return true;
	}
	else return false;
}

bool DoubleKeyHashTable::wennDiagramTwoLists(vector<UINT64>* list1, vector<UINT64>* list2, vector<UINT64>* list1only, vector<UINT64>* list2only, vector<UINT64>* list12)//sorted from smaller to larger ID in the list
{
	unsigned int i,j;
	i=0;j=0;
	while(list1!=NULL&&list2!=NULL&& i<list1->size()&&j<list2->size())
	{
		if(list1->at(i)<list2->at(j))
		{
			UINT64 smallvalue = list1->at(i);
			list1only->push_back(smallvalue);
			i++;
		}
		else if(list1->at(i)>list2->at(j))
		{
			UINT64 smallvalue = list2->at(j);
			list2only->push_back(smallvalue);
			j++;
		}
		else //same value in both list
		{
			UINT64 value = list1->at(i);
			list12->push_back(value);
			i++;j++;
		}
	}

	while(list1!=NULL && i<list1->size())
	{
		UINT64 smallvalue = list1->at(i);
		list1only->push_back(smallvalue);
		i++;
	}

	while(list2!=NULL && j<list2->size())
	{
		UINT64 smallvalue = list2->at(j);
		list2only->push_back(smallvalue);
		j++;
	}
	return true;
}

//keymatchmode = 1(only left key matched),2(only right key matched),3(both matched)
bool DoubleKeyHashTable::createAlignment(Alignment* subjectAlignment, UINT8 querymode , UINT8 keymatchmode, UINT64 subjectKeyStart)
{
	QueryRead* queryRead = subjectAlignment->queryRead;
	SubjectRead* subjectRead = subjectAlignment->subjectRead;
	string subjectString = subjectRead->getSequence();
	string queryString;
	int hashKeyLength = this->hashKeyLength_left+this->hashKeyLength_right;
	int start,end =0;
	int seedStart, seedEnd = 0;
	switch (querymode) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
	{
	case 0:
		subjectAlignment->queryOrientation = true;
		subjectAlignment->subjectStart = -subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		queryString = queryRead->getSequence();
		start = 0;
		end = (subjectAlignment->subjectEnd <= subjectAlignment->queryEnd)?subjectAlignment->subjectEnd:subjectAlignment->queryEnd;
		switch(keymatchmode)
		{
		case 1:
			seedStart = 0;
			seedEnd = this->hashKeyLength_left-1;
			break;
		case 2:
			seedStart = this->hashKeyLength_left;
			seedEnd = hashKeyLength-1;
			break;
		case 3:
			seedStart = 0;
			seedEnd = hashKeyLength-1;
			break;
		default: return false;
		}

		break;
	case 1:
		subjectAlignment->queryOrientation = true;
		subjectAlignment->subjectStart = queryRead->getReadLength()-hashKeyLength-subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		queryString = queryRead->getSequence();
		start = (0>=subjectAlignment->subjectStart)?0:subjectAlignment->subjectStart;
		end = subjectAlignment->queryEnd;
		switch(keymatchmode)
		{
		case 1:
			seedStart = subjectAlignment->queryEnd-hashKeyLength+1;
			seedEnd = subjectAlignment->queryEnd-this->hashKeyLength_right;
			break;
		case 2:
			seedStart = subjectAlignment->queryEnd-this->hashKeyLength_right+1;
			seedEnd = subjectAlignment->queryEnd;
			break;
		case 3:
			seedStart = subjectAlignment->queryEnd-hashKeyLength+1;
			seedEnd = subjectAlignment->queryEnd;
			break;
		default: return false;
		}
		break;
	case 2:
		subjectAlignment->queryOrientation = false;
		subjectAlignment->subjectStart = -subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		queryString = queryRead->reverseComplement();
		start = 0;
		end = (subjectAlignment->subjectEnd <= subjectAlignment->queryEnd)?subjectAlignment->subjectEnd:subjectAlignment->queryEnd;
		switch(keymatchmode)
		{
		case 1:
			seedStart = 0;
			seedEnd = this->hashKeyLength_left-1;
			break;
		case 2:
			seedStart = this->hashKeyLength_left;
			seedEnd = hashKeyLength-1;
			break;
		case 3:
			seedStart = 0;
			seedEnd = hashKeyLength-1;
			break;
		default: return false;
		}
		break;
	case 3:
		subjectAlignment->queryOrientation = false;
		subjectAlignment->subjectStart = queryRead->getReadLength()-hashKeyLength-subjectKeyStart;
		subjectAlignment->queryEnd = queryRead->getReadLength()-1;
		subjectAlignment->subjectEnd = subjectAlignment->subjectStart + subjectRead->getReadLength() -1;
		queryString = queryRead->reverseComplement();
		start = (0>=subjectAlignment->subjectStart)?0:subjectAlignment->subjectStart;
		end = subjectAlignment->queryEnd;
		switch(keymatchmode)
		{
		case 1:
			seedStart = subjectAlignment->queryEnd-hashKeyLength+1;
			seedEnd = subjectAlignment->queryEnd-this->hashKeyLength_right;
			break;
		case 2:
			seedStart = subjectAlignment->queryEnd-this->hashKeyLength_right+1;
			seedEnd = subjectAlignment->queryEnd;
			break;
		case 3:
			seedStart = subjectAlignment->queryEnd-hashKeyLength+1;
			seedEnd = subjectAlignment->queryEnd;
			break;
		default: return false;
		}
		break;
	default: return false;
	}
	bool alignsuccess = this->doAlignmentWithSeed(subjectAlignment,queryString, subjectString,start, end, seedStart, seedEnd);
	return alignsuccess;
}

bool DoubleKeyHashTable::searchHashTable(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subLeftString, subRightString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()- this->hashKeyLength_left-this->hashKeyLength_right; j++)
	{
		subLeftString = subjectReadString.substr(j,this->hashKeyLength_left);
		subRightString = subjectReadString.substr(j+this->hashKeyLength_left, this->hashKeyLength_right);
		vector<UINT64> * listOfReads_left=this->getListOfReads(subLeftString);
		vector<UINT64> * listOfReads_right=this->getListOfReads(subRightString);
		vector<UINT64> * listOfReadsList[8];
		for(int i=0;i<8;i++)
			listOfReadsList[i] = new vector<UINT64>();

		if(listOfReads_left!=NULL && !listOfReads_left->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_left->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_left->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 0 || queryMode == 1 || queryMode == 2 || queryMode == 3)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()>queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}
		if(listOfReads_right!=NULL && !listOfReads_right->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_right->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_right->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 4 || queryMode == 5 || queryMode == 6 || queryMode == 7)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()>queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}

		for(UINT8 i =0;i<=3;i++)
		{
			if(listOfReadsList[i]->empty()==false || listOfReadsList[i+4]->empty()==false)
			{
			vector<UINT64>* KeyLeftOnlyList = new vector<UINT64>;
			vector<UINT64>* KeyRightOnlyList = new vector<UINT64>;
			vector<UINT64>* BothKeyList = new vector<UINT64>;
			this->wennDiagramTwoLists(listOfReadsList[i],listOfReadsList[i+4],KeyLeftOnlyList, KeyRightOnlyList, BothKeyList);

				for(unsigned int k=0;k<BothKeyList->size();k++)
				{
					UINT64 currentID = BothKeyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 3;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;

				}

				for(unsigned int k=0;k<KeyLeftOnlyList->size();k++)
				{
					UINT64 currentID = KeyLeftOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 1;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
				for(unsigned int k=0;k<KeyRightOnlyList->size();k++)
				{
					UINT64 currentID = KeyRightOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 2;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
				KeyLeftOnlyList->clear();
				KeyRightOnlyList->clear();
				BothKeyList->clear();
				delete KeyLeftOnlyList;
				delete KeyRightOnlyList;
				delete BothKeyList;
			}// end if


		}//end for

		for(int i=0;i<8;i++)
		{
			listOfReadsList[i]->clear();
			delete listOfReadsList[i];
		}
	}
	return true;
}

bool DoubleKeyHashTable::searchHashTableFullScan(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subLeftString, subRightString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength_left-this->hashKeyLength_right; j++)
	{
		subLeftString = subjectReadString.substr(j,this->hashKeyLength_left);
		subRightString = subjectReadString.substr(j+this->hashKeyLength_left, this->hashKeyLength_right);
		vector<UINT64> * listOfReads_left=this->getListOfReads(subLeftString);
		vector<UINT64> * listOfReads_right=this->getListOfReads(subRightString);
		vector<UINT64> * listOfReadsList[8];
		for(int i=0;i<8;i++)
			listOfReadsList[i] = new vector<UINT64>();

		if(listOfReads_left!=NULL && !listOfReads_left->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_left->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_left->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 0 || queryMode == 1 || queryMode == 2 || queryMode == 3)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()!=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}
		if(listOfReads_right!=NULL && !listOfReads_right->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_right->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_right->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 4 || queryMode == 5 || queryMode == 6 || queryMode == 7)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()!=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}

		for(UINT8 i =0;i<=3;i++)
		{
			if(listOfReadsList[i]->empty()==false || listOfReadsList[i+4]->empty()==false)
			{
			vector<UINT64>* KeyLeftOnlyList = new vector<UINT64>;
			vector<UINT64>* KeyRightOnlyList = new vector<UINT64>;
			vector<UINT64>* BothKeyList = new vector<UINT64>;
			this->wennDiagramTwoLists(listOfReadsList[i],listOfReadsList[i+4],KeyLeftOnlyList, KeyRightOnlyList, BothKeyList);

				for(unsigned int k=0;k<BothKeyList->size();k++)
				{
					UINT64 currentID = BothKeyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 3;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;

				}

				for(unsigned int k=0;k<KeyLeftOnlyList->size();k++)
				{
					UINT64 currentID = KeyLeftOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 1;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
				for(unsigned int k=0;k<KeyRightOnlyList->size();k++)
				{
					UINT64 currentID = KeyRightOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 2;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
				KeyLeftOnlyList->clear();
				KeyRightOnlyList->clear();
				BothKeyList->clear();
				delete KeyLeftOnlyList;
				delete KeyRightOnlyList;
				delete BothKeyList;
			}// end if


		}//end for

		for(int i=0;i<8;i++)
		{
			listOfReadsList[i]->clear();
			delete listOfReadsList[i];
		}
	}
	return true;
}

bool DoubleKeyHashTable::searchHashTableFullScanAndMarkQueryRead(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subLeftString, subRightString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength_left-this->hashKeyLength_right; j++)
	{
		subLeftString = subjectReadString.substr(j,this->hashKeyLength_left);
		subRightString = subjectReadString.substr(j+this->hashKeyLength_left, this->hashKeyLength_right);
		vector<UINT64> * listOfReads_left=this->getListOfReads(subLeftString);
		vector<UINT64> * listOfReads_right=this->getListOfReads(subRightString);
		vector<UINT64> * listOfReadsList[8];
		for(int i=0;i<8;i++)
			listOfReadsList[i] = new vector<UINT64>();

		if(listOfReads_left!=NULL && !listOfReads_left->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_left->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_left->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 0 || queryMode == 1 || queryMode == 2 || queryMode == 3)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()!=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}
		if(listOfReads_right!=NULL && !listOfReads_right->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_right->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_right->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 4 || queryMode == 5 || queryMode == 6 || queryMode == 7)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()!=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}

		for(UINT8 i =0;i<=3;i++)
		{
			if(listOfReadsList[i]->empty()==false || listOfReadsList[i+4]->empty()==false)
			{
			vector<UINT64>* KeyLeftOnlyList = new vector<UINT64>;
			vector<UINT64>* KeyRightOnlyList = new vector<UINT64>;
			vector<UINT64>* BothKeyList = new vector<UINT64>;
			this->wennDiagramTwoLists(listOfReadsList[i],listOfReadsList[i+4],KeyLeftOnlyList, KeyRightOnlyList, BothKeyList);

				for(unsigned int k=0;k<BothKeyList->size();k++)
				{
					UINT64 currentID = BothKeyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();

					if(queryRead->flag4Removal==false)
					{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 3;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
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

				for(unsigned int k=0;k<KeyLeftOnlyList->size();k++)
				{
					UINT64 currentID = KeyLeftOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();

					if(queryRead->flag4Removal==false)
					{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 1;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
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
				for(unsigned int k=0;k<KeyRightOnlyList->size();k++)
				{
					UINT64 currentID = KeyRightOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();

					if(queryRead->flag4Removal==false)
					{
					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 2;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
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
				KeyLeftOnlyList->clear();
				KeyRightOnlyList->clear();
				BothKeyList->clear();
				delete KeyLeftOnlyList;
				delete KeyRightOnlyList;
				delete BothKeyList;
			}// end if


		}//end for

		for(int i=0;i<8;i++)
		{
			listOfReadsList[i]->clear();
			delete listOfReadsList[i];
		}
	}
	return true;
}

bool DoubleKeyHashTable::searchHashTable(SubjectEdge * subjectEdge, bool fullScan)
{
	if(fullScan == false) return this->searchHashTable(subjectEdge);
	else return this->searchHashTableFullScan(subjectEdge);
}

bool DoubleKeyHashTable::searchHashTableFullScanAndRecordKeys(SubjectEdge * subjectEdge)
{
	SubjectRead *subjectRead = subjectEdge->subjectRead; 	// Get the current read subject read.
	string subjectReadString = subjectRead->getSequence(); 		// Get the forward string of subject read.
	string subLeftString, subRightString;
	for(INT64 j = 0; j <= (INT64)subjectRead->getReadLength()-this->hashKeyLength_left-this->hashKeyLength_right; j++)
	{
		subLeftString = subjectReadString.substr(j,this->hashKeyLength_left);
		subRightString = subjectReadString.substr(j+this->hashKeyLength_left, this->hashKeyLength_right);

		UINT64 hashValue;
		vector<UINT64> * listOfReads_left=this->getListOfReads(subLeftString, hashValue); // Search the string in the hash table.
		subjectEdge->addKeyToMap(hashValue,j); //only need to record the left string hashvalue.
		vector<UINT64> * listOfReads_right=this->getListOfReads(subRightString);
		vector<UINT64> * listOfReadsList[8];
		for(int i=0;i<8;i++)
			listOfReadsList[i] = new vector<UINT64>();

		if(listOfReads_left!=NULL && !listOfReads_left->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_left->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_left->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 0 || queryMode == 1 || queryMode == 2 || queryMode == 3)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()!=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}
		if(listOfReads_right!=NULL && !listOfReads_right->empty())
		{
			for(INT64 k = 0; k < (INT64)listOfReads_right->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads_right->at(k);
				UINT64 queryReadID = data & 0X1FFFFFFFFFFFFFFF;
				UINT64 queryMode = data >> this->numberOfLSB;
				if(queryMode == 4 || queryMode == 5 || queryMode == 6 || queryMode == 7)
				{
					QueryRead *queryRead = this->dataSet->getReadFromID(queryReadID);
					if(subjectRead->getName()!=queryRead->getName()) 			// No need to discover the same edge again. Only need to explore half of the combinations
					{
						listOfReadsList[queryMode]->push_back(queryReadID);
					}
				}
			}
		}

		for(UINT8 i =0;i<=3;i++)
		{
			if(listOfReadsList[i]->empty()==false || listOfReadsList[i+4]->empty()==false)
			{
			vector<UINT64>* KeyLeftOnlyList = new vector<UINT64>;
			vector<UINT64>* KeyRightOnlyList = new vector<UINT64>;
			vector<UINT64>* BothKeyList = new vector<UINT64>;
			this->wennDiagramTwoLists(listOfReadsList[i],listOfReadsList[i+4],KeyLeftOnlyList, KeyRightOnlyList, BothKeyList);

				for(unsigned int k=0;k<BothKeyList->size();k++)
				{
					UINT64 currentID = BothKeyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 3;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;

				}

				for(unsigned int k=0;k<KeyLeftOnlyList->size();k++)
				{
					UINT64 currentID = KeyLeftOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 1;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
				for(unsigned int k=0;k<KeyRightOnlyList->size();k++)
				{
					UINT64 currentID = KeyRightOnlyList->at(k);
					QueryRead* queryRead = this->dataSet->getReadFromID(currentID);
//					string queryString = (i==0 || i==1)?queryRead->getSequence():queryRead->reverseComplement();


					Alignment *subjectAlignment = new Alignment(subjectRead, queryRead);
					UINT8 keyMatchMode = 2;

					if(this->createAlignment(subjectAlignment, i, keyMatchMode, j))
					{
					subjectEdge->addAlignment(subjectAlignment);
					}
					else
						delete subjectAlignment;
				}
				KeyLeftOnlyList->clear();
				KeyRightOnlyList->clear();
				BothKeyList->clear();
				delete KeyLeftOnlyList;
				delete KeyRightOnlyList;
				delete BothKeyList;
			}// end if


		}//end for

		for(int i=0;i<8;i++)
		{
			listOfReadsList[i]->clear();
			delete listOfReadsList[i];
		}
	}
	return true;
}

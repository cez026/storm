/*
 * HashTableMethod.h
 *
 *  Created on: Feb 13, 2015
 *      Author: qy2
 */

#ifndef HASHTABLEMETHOD_H_
#define HASHTABLEMETHOD_H_

#include "Config.h"
#include "QueryDataset.h"
#include "HashTable.h"
#include "QueryRead.h"
#include "Alignment.h"

class HashTableMethod {
protected:
	QueryDataset *dataSet;
	HashTable* hashTable;
	UINT64 numberOfHashCollision;				// Counted total number of hash collisions. For debugging only.
	UINT64 maxSingleHashCollision;				// Counted maximal number of hash collision for a single case.

public:
	HashTableMethod();
	virtual ~HashTableMethod();

	virtual bool createHashTables()=0;
	virtual string getReadSubstring(UINT64 readID, UINT8 mode)=0;
	virtual bool insertQueryDataset(QueryDataset* d)=0;
	virtual bool insertQueryDataset_leftReadForMerge(QueryDataset* d)=0;
	virtual bool insertQueryDataset_rightAlignOnly(QueryDataset* d)=0;
	virtual bool insertQueryRead(QueryRead *read,string subString, UINT8 mode)=0;
	virtual vector<UINT64> * getListOfReads(string subString)=0;
	virtual vector<UINT64> * getListOfReads(string subString, UINT64& hashValue)=0;
	virtual bool searchHashTable(SubjectEdge * subjectEdge)=0;
	virtual bool searchHashTableFullScan(SubjectEdge * subjectEdge)=0;
	virtual bool searchHashTable(SubjectEdge * subjectEdge, bool fullScan)=0;
	virtual bool searchHashTableFullScanAndMarkQueryRead(SubjectEdge * subjectEdge)=0;
	virtual bool searchHashTableFullScanAndRecordKeys(SubjectEdge * subjectEdge)=0;
	QueryDataset* getDataset(){return this->dataSet;};
	HashTable* getHashTable(){return this->hashTable;};
};

#endif /* HASHTABLEMETHOD_H_ */

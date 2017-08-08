/*
 * HashTable.h
 *
 *  Created on: Jan 6, 2015
 *      Author: qy2
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "Config.h"
#include "QueryDataset.h"


//Since there is no internal keystring stored, we cannot handle the collision internally
class HashTable {

private:

	UINT64 hashTableSize; 					// Size of hash table, primer number usually.



	vector < vector<UINT64> *> *hashTable; 		// Main structure of hashtable, storing read identifiers.


	UINT64 getPrimeLargerThanNumber(UINT64 number);
	void setHashTableSizeAndInitialize(UINT64 size);


public:

	HashTable(UINT64 size);
	~HashTable();
	UINT64 hashFunction(const string & subString);

	UINT64 getHashTableSize();


	bool isEmptyAt(UINT64 hashTableIndex);
	vector<UINT64> * getReadIDListAt(UINT64 hashTableIndex);
	bool insertReadIDAt(UINT64 hashTableIndex, UINT64 readID);



};

#endif /* HASHTABLE_H_ */

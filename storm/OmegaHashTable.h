/*
 * HashTable.h
 *
 *  Created on: Apr 22, 2013
 *      Author: b72
 */

#ifndef OMEGAHASHTABLE_H_
#define OMEGAHASHTABLE_H_
#include "Config.h"
#include "QueryDataset.h"

/**********************************************************************************************************************
	Class to store hashtable.
**********************************************************************************************************************/
class OmegaHashTable{
	private:
		QueryDataset *dataSet;							// Pointer of the dataset.
		UINT64 hashTableSize; 						// Ted: Size of the hash table. This is the prime number of mod operation.
		vector < vector<UINT64> *> *hashTable; 		// Ted: List of hash keys of read number (ID) and orientation.
		UINT16 hashStringLength;					// Ted: Length of prefix and suffix of the reads to hash. This is equal to the minumum overlap length.
		UINT64 numberOfHashCollision;				// Counter to count the number of hash collisions. For debugging only.
		bool insertIntoTable(QueryRead *read, string substring, UINT64 orientation);	// Insert a string in the hash table.
		bool hashRead(QueryRead *read); 					// Ted: Hash prefix and suffix of the read and its reverse complement in the hash table. Turn over to the constant
		bool hashRead_half(QueryRead *read); 					//Hash prefix and suffix of the read and its reverse complement in the hash table. Turn over to the constant
		bool hashRead_left(QueryRead *read); //only hash the left paired end read, with suffix of the forward, and prefix of the reverse
		void setHashTableSize(UINT64 size); 		// Set the size of the hash table.
	public:
		OmegaHashTable(void);							// Default constructor.
		OmegaHashTable(QueryDataset *d);						// Another constructor.
		~OmegaHashTable();								// Destructor.
		bool insertDataset(QueryDataset *d, UINT64 minOverlapLength);	// Insert the dataset in the hash table.
		bool insertDataset_half(QueryDataset *d, UINT64 minOverlapLength);	// Insert the dataset in the hash table.
		bool insertDataset_left(QueryDataset * d, UINT64 minOverlapLength); // Insert the dataset (only left part of the paired end reads) in the hash table.
		bool insertDataset_rightAlignOnly(QueryDataset* d, UINT64 minOverlapLength);
		vector<UINT64> * getListOfReads(string subString); 			// Get the list of reads that contain subString as prefix or suffix.
		UINT64 hashFunction(string subString); 						// Hash function.
		UINT64 getHashTableSize(void){return hashTableSize;}		// Get the size of the hash table.
		UINT64 getHashStringLength() {return hashStringLength;}		// Get the hash string length.
		QueryDataset * getDataset(void){return dataSet;}					// Get the pointer to the dataset.
};


#endif /* OMEGAHASHTABLE_H_ */

/*
 * QueryDataset.h
 *
 *  Created on: Dec 23, 2014
 *      Author: qy2
 */

#ifndef QUERYDATASET_H_
#define QUERYDATASET_H_

#include "Config.h"
#include "QueryRead.h"
#include "PairedEndRead.h"

class QueryDataset {

	UINT64 numberOfReads;								// Number of total reads present in the dataset.
	UINT64 numberOfUniqueReads; 						// number of unique reads in the dataset.

	UINT64 shortestReadLength;
	UINT64 longestReadLength;

	UINT16 dataset_minimumoverlaplength;
	string dataset_QueryFilename;
	UINT16 dataset_leftminimumoverlaplength;
	UINT16 dataset_rightminimumoverlaplength;

	vector<QueryRead *>* queryReadList;

	vector<PairedEndRead*>* pairedEndReadList;


	bool buildDatasetFromMatePairFileFF(const string& QueryFilename);
	bool buildDatasetFromMatePairFileFR(const string& QueryFilename);
	bool buildDatasetFromMatePairFileRF(const string& QueryFilename);
	bool buildDatasetFromMatePairFileRR(const string& QueryFilename);

	bool duplicateFilter();//filter the identical reads and add the number to frequency
	void sortReads();//in order to assign IDs and facilitate binary sequence search


public:
	QueryDataset();
	QueryDataset(const string & QueryFilename);
	~QueryDataset();


	static bool qualityFilter(string & sequence);

	bool buildDataset(const string & QueryFilename);
	bool buildDatasetParallel(const string & QueryFilename, UINT16 numberOfThreads);

	bool buildDatasetNoSortingDuplicateRemoving(const string & QueryFilename, UINT16 numberOfThreads);

	bool buildDatasetFromMatePairFile(const string& QueryFilename, UINT8 orient);
	bool buildDatasetNoSortingDuplicateRemoving();
	bool buildDataset();
	UINT64 getNumberOfReads(); 						// Get the number of total reads in the database.
	UINT64 getNumberOfUniqueReads(); 				// Get the number of unique reads in the database.
	vector<PairedEndRead *>* getPairedEndReadList();

//	QueryRead * getReadFromString(const string & read); 		// Find a read in the database given the string. Uses binary search in the list of reads.
	QueryRead * getReadFromID(UINT64 ID); 					// Find a read in the database given the ID in constant time.

	//for debugging, not implemented yet
//	bool printDataset(int from, int to); 					// Print few the reads in the dataset. For debuggin only.
//	void saveReads(string fileName); // Save all the sorted unique reads in a text file. Used for debugging.

};

#endif /* QUERYDATASET_H_ */

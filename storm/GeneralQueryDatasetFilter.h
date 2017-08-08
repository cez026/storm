/*
 * GeneralQueryDatasetFilter.h
 *
 *  Created on: Apr 29, 2015
 *      Author: qy2
 */

#ifndef GENERALQUERYDATASETFILTER_H_
#define GENERALQUERYDATASETFILTER_H_
#include "Config.h"
#include "SingleKeyHashTable.h"
#include "DoubleKeyHashTable.h"
#include "QueryDataset.h"
#include "HashTableMethod.h"
#include "SubjectDataset.h"
#include "Alignment.h"
class GeneralQueryDatasetFilter {

	ofstream filePointer;
	HashTableMethod * hashTableMethod;

	UINT16 minimumOverlapLength;


public:
	GeneralQueryDatasetFilter(HashTableMethod * hashTableMethod);
	~GeneralQueryDatasetFilter();

	bool start();
	bool searchHashTableAndMarkQueryReads(SubjectEdge * subjectEdge);
	bool checkOverlapForContainedRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start);
	bool checkIdenticalRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start);
	void printToFile();
};

#endif /* GENERALQUERYDATASETFILTER_H_ */

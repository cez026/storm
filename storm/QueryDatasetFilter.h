/*
 * QueryDatasetFilter.h
 *
 *  Created on: Feb 23, 2015
 *      Author: qy2
 */

#ifndef QUERYDATASETFILTER_H_
#define QUERYDATASETFILTER_H_

#include "Config.h"
#include "OmegaHashTable.h"
#include "SubjectDataset.h"

class QueryDatasetFilter {
	OmegaHashTable * omegaHashTable;
	QueryDataset *dataset;


	//tag				tag			frequency
	//unique reads		0			1
	//with max name		0			>1
	//not max name		1			>1
	//contained reads	+2
//	vector<UINT8> * tagList;
	//store the super reads name of this contained read
//	vector<string> * superNameList;
	//store the length of the super containing reads
//	vector<int> * superReadLength;

	ofstream filePointer;
public:
	QueryDatasetFilter(OmegaHashTable * omegaHashTable);
	~QueryDatasetFilter();
	bool start();
	bool start_backup();
	bool searchHashTable(SubjectEdge * subjectEdge);
	bool checkOverlapForContainedRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start);
	bool checkIdenticalRead(SubjectRead *read1, QueryRead *read2, UINT64 orient, UINT64 start);
	void printToFile();
};

#endif /* QUERYDATASETFILTER_H_ */
